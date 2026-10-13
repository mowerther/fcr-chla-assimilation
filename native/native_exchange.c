/* SPDX-License-Identifier: GPL-3.0-or-later
 * State exchange for the published GLM-AED executable (Falling Creek Reservoir).
 * Loaded with LD_PRELOAD. Interposes the write of the hourly time coordinate and,
 * at each scheduled event hour, writes the forecast state, reads the updated
 * temperature and phytoplankton, writes the accepted state and waits for a
 * continue signal. No integration routine is replaced.
 * Layout and signatures follow GLM3 commit 8453b1d (glm.h, glm_types.h, GPL-3.0);
 * runtime validation pins this layout.
 */
#define _GNU_SOURCE
#include <dlfcn.h>
#include <netcdf.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <sys/inotify.h>
#include <poll.h>
#include "glm.h"
#include "glm_types.h"

static void die(const char *s) { fprintf(stderr,"exchange: %s\n",s); exit(91); }
static void *sym(const char *s) { void *p=dlsym(RTLD_DEFAULT,s); if(!p) die(s); return p; }
static void wait_file(const char *p) {
    /* Event-driven wait, so that paused members do not poll.
     * The watch is installed before the second existence query, so no wake is lost.
     */
    static int fd=-1;
    if(access(p,F_OK)==0)return;
    if(fd<0) {
        char dir[1024];snprintf(dir,sizeof(dir),"%s",p);
        char *slash=strrchr(dir,'/');if(!slash)die("exchange directory");*slash=0;
        fd=inotify_init1(IN_CLOEXEC);
        if(fd<0||inotify_add_watch(fd,dir,IN_CREATE|IN_MOVED_TO)<0)die("inotify setup");
    }
    while(access(p,F_OK)) {
        struct pollfd f={fd,POLLIN,0};
        if(poll(&f,1,3600000)<=0)die("exchange timeout");
        char events[8192];if(read(fd,events,sizeof(events))<=0)die("exchange notification");
    }
}
static void path(char *p,const char *root,const char *kind,int rec) { snprintf(p,1024,"%s/%s_%03d.bin",root,kind,rec); }

static void dump(const char *root,const char *kind,int rec,int n,int nv,int ml,LakeDataType *lake,double *wq,int *ids) {
    char p[1024],tmp[1060]; path(p,root,kind,rec); snprintf(tmp,sizeof(tmp),"%s.partial",p);
    FILE *f=fopen(tmp,"wb"); if(!f)die("snapshot open");
    int32_t header[8]={0x46435258,rec,n,nv,ml,ids[0],ids[1],ids[2]};
    fwrite(header,4,8,f);
    for(int i=0;i<n;i++) {
        double row[7]={lake[i].Height,lake[i].LayerVol,lake[i].Temp,lake[i].Salinity,lake[i].Density,lake[i].LayerArea,lake[i].ExtcCoefSW};
        fwrite(row,8,7,f);
    }
    for(int v=0;v<nv;v++) fwrite(wq+v*ml,8,n,f);
    if(fclose(f))die("snapshot close");
    if(rename(tmp,p))die("snapshot rename");
}

static void exchange(int rec) {
    const char *root=getenv("GLM_EXCHANGE_DIR"); if(!root)return;
    int n=*(int*)sym("NumLayers"),nv=*(int*)sym("Num_WQ_Vars"),ml=*(int*)sym("MaxLayers");
    LakeDataType *lake=*(LakeDataType**)sym("Lake"); double *wq=*(double**)sym("WQ_Vars");
    int (*indexfn)(const char*,size_t*)=sym("aed_var_index_c");
    double (*density)(double,double)=sym("calculate_density");
    const char *names[]={"PHY_cyano","PHY_green","PHY_diatom"}; int ids[3];
    for(int j=0;j<3;j++){size_t len=strlen(names[j]);ids[j]=indexfn(names[j],&len);if(ids[j]<0||ids[j]>=nv)die("phyto index");}
    if(n<2||n>ml||ml!=500||nv!=24||lake[n-1].Height<8||lake[n-1].Height>11)die("ABI validation");
    dump(root,"forecast",rec,n,nv,ml,lake,wq,ids);
    if(getenv("GLM_CAPTURE_ONLY"))return;
    char p[1024];path(p,root,"update",rec);wait_file(p);
    FILE *f=fopen(p,"rb"); int32_t header[2];if(fread(header,4,2,f)!=2||header[0]!=rec||header[1]!=n)die("update header");
    double *a=malloc(n*4*sizeof(double));if(fread(a,8,n*4,f)!=(size_t)n*4)die("update length");fclose(f);
    for(int i=0;i<n;i++) {
        if(!isfinite(a[i*4])||a[i*4]<-2||a[i*4]>45)die("temperature bound");
        for(int j=1;j<4;j++)if(!isfinite(a[i*4+j])||a[i*4+j]<0.03-1e-10||a[i*4+j]>10000)die("phyto bound");
    }
    for(int i=0;i<n;i++) {
        /* Density is recalculated only where temperature changed, so unchanged memory stays exact. */
        if(a[i*4]!=lake[i].Temp){lake[i].Temp=a[i*4];lake[i].Density=density(lake[i].Temp,lake[i].Salinity);}
        for(int j=0;j<3;j++)wq[ids[j]*ml+i]=a[i*4+j+1];
    }
    free(a);dump(root,"accepted",rec,n,nv,ml,lake,wq,ids);
    path(p,root,"continue",rec);wait_file(p);
}

int nc_put_vara(int ncid,int varid,const size_t *start,const size_t *count,const void *values) {
    static int (*original)(int,int,const size_t*,const size_t*,const void*)=NULL;
    static int timeid=-1,activeid=-1,last_hour=-1;
    static int hours[2000], nhours=-1, next=0;
    if(!original){original=dlsym(RTLD_NEXT,"nc_put_vara");if(!original)die("nc_put_vara resolution");}
    if(getenv("GLM_EXCHANGE_DIR")) {
        if(nhours<0) {
            const char *schedule=getenv("GLM_EVENT_HOURS");if(!schedule)die("missing hour schedule");
            FILE *f=fopen(schedule,"r");if(!f)die("schedule open");nhours=0;
            while(nhours<2000 && fscanf(f,"%d",&hours[nhours])==1)nhours++;
            fclose(f);if(nhours==0||nhours==2000)die("schedule length");
            for(int i=1;i<nhours;i++)if(hours[i]<=hours[i-1])die("schedule order");
        }
        if(activeid!=ncid||timeid<0){activeid=ncid;if(nc_inq_varid(ncid,"time",&timeid)!=NC_NOERR)timeid=-1;}
        if(varid==timeid) {
            double h=*(const double*)values;int ih=(int)llround(h);
            if(fabs(h-ih)>1e-8)die("hourly native clock expected");
            if(ih!=last_hour) {
                last_hour=ih;
                if(next<nhours && ih>hours[next])die("native schedule skipped");
                if(next<nhours && ih==hours[next]) { exchange(next);next++; }
            }
        }
    }
    /* With GLM_SKIP_OUTPUT the netCDF output is not written. The accepted
     * native states in the exchange files are the outputs of a run.
     */
    if(getenv("GLM_SKIP_OUTPUT"))return NC_NOERR;
    /* With GLM_COMPACT_OUTPUT only the latest output record is kept. */
    size_t compact_start[NC_MAX_VAR_DIMS]; const size_t *outstart=start;
    if(getenv("GLM_COMPACT_OUTPUT")) {
        int nd,dims[NC_MAX_VAR_DIMS];char name[NC_MAX_NAME+1];
        if(nc_inq_varndims(ncid,varid,&nd)==NC_NOERR && nd>0 &&
           nc_inq_vardimid(ncid,varid,dims)==NC_NOERR &&
           nc_inq_dimname(ncid,dims[0],name)==NC_NOERR && !strcmp(name,"time")) {
            memcpy(compact_start,start,nd*sizeof(size_t));compact_start[0]=0;outstart=compact_start;
        }
    }
    int result=original(ncid,varid,outstart,count,values);
    return result;
}
