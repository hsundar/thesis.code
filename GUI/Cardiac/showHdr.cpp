#include #include "dbh.h"

void ShowHdr(char *, struct dsr *); 
void swap_long(unsigned char *); 
void swap_short(unsigned char *);

main(argc,argv) int argc; char **argv;
{ 
  struct dsr hdr; 
  int size; 
  double cmax, 
         cmin; FILE *fp;

  if((fp=fopen(argv[1],"r"))==NULL)
  {
    fprintf(stderr,"Can't open:\n", argv[1]); exit(0);
  } fread(&hdr,1,sizeof(struct dsr),fp);

  if(hdr.dime.dim[0] 15)
    swap_hdr(&hdr);

  ShowHdr(argv[1], &hdr);


}




void ShowHdr(fileName,hdr) struct dsr *hdr; char *fileName; { int i; char
  string[128]; printf("Analyze Header Dump of: \n", fileName); /* Header Key */
  printf("sizeof_hdr: \n", hdr->hk.sizeof_hdr); printf("data_type: \n",
      hdr->hk.data_type); printf("db_name: \n", hdr->hk.db_name); printf("extents:
        \n", hdr->hk.extents); printf("session_error: \n", hdr->hk.session_error);
  printf("regular: \n", hdr->hk.regular); printf("hkey_un0: \n",
      hdr->hk.hkey_un0);

  /* Image Dimension */ for(i=0;i<8;i++)
    printf("dim[%d]: \n", i, hdr->dime.dim[i]);

  strncpy(string,hdr->dime.vox_units,4); printf("vox_units: \n", string);

  strncpy(string,hdr->dime.cal_units,8); printf("cal_units: \n", string);
  printf("unused1: \n", hdr->dime.unused1); printf("datatype: \n",
      hdr->dime.datatype); printf("bitpix: \n", hdr->dime.bitpix);

  for(i=0;i<8;i++)
    printf("pixdim[%d]: \n",i, hdr->dime.pixdim[i]);

  printf("vox_offset: \n", hdr->dime.vox_offset); printf("funused1: \n",
      hdr->dime.funused1); printf("funused2: \n", hdr->dime.funused2);
  printf("funused3: \n", hdr->dime.funused3); printf("cal_max: \n",
      hdr->dime.cal_max); printf("cal_min: \n", hdr->dime.cal_min);
  printf("compressed: \n", hdr->dime.compressed); printf("verified: \n",
      hdr->dime.verified); printf("glmax: \n", hdr->dime.glmax); printf("glmin: \n",
        hdr->dime.glmin);

  /* Data History */ strncpy(string,hdr->hist.descrip,80); printf("descrip: \n",
      string); strncpy(string,hdr->hist.aux_file,24); printf("aux_file: \n", string);
  printf("orient: \n", hdr->hist.orient);

  strncpy(string,hdr->hist.originator,10); printf("originator: \n", string);

  strncpy(string,hdr->hist.generated,10); printf("generated: \n", string);


  strncpy(string,hdr->hist.scannum,10); printf("scannum: \n", string);

  strncpy(string,hdr->hist.patient_id,10); printf("patient_id: \n", string);

  strncpy(string,hdr->hist.exp_date,10); printf("exp_date: \n", string);

  strncpy(string,hdr->hist.exp_time,10); printf("exp_time: \n", string);

  strncpy(string,hdr->hist.hist_un0,10); printf("hist_un0: \n", string);

  printf("views: \n", hdr->hist.views); printf("vols_added: \n",
      hdr->hist.vols_added); printf("start_field: \n", hdr->hist.start_field);
  printf("field_skip: \n", hdr->hist.field_skip); printf("omax: \n",
      hdr->hist.omax); printf("omin: \n", hdr->hist.omin); printf("smin: \n",
        hdr->hist.smax); printf("smin: \n", hdr->hist.smin);

}


swap_hdr(pntr) struct dsr *pntr;
{ swap_long(&pntr->hk.sizeof_hdr) ; swap_long(&pntr->hk.extents) ;
  swap_short(&pntr->hk.session_error) ; swap_short(&pntr->dime.dim[0]) ;
  swap_short(&pntr->dime.dim[1]) ; swap_short(&pntr->dime.dim[2]) ;
  swap_short(&pntr->dime.dim[3]) ; swap_short(&pntr->dime.dim[4]) ;
  swap_short(&pntr->dime.dim[5]) ; swap_short(&pntr->dime.dim[6]) ;
  swap_short(&pntr->dime.dim[7]) ; swap_short(&pntr->dime.unused1) ;
  swap_short(&pntr->dime.datatype) ; swap_short(&pntr->dime.bitpix) ;
  swap_long(&pntr->dime.pixdim[0]) ; swap_long(&pntr->dime.pixdim[1]) ;
  swap_long(&pntr->dime.pixdim[2]) ; swap_long(&pntr->dime.pixdim[3]) ;
  swap_long(&pntr->dime.pixdim[4]) ; swap_long(&pntr->dime.pixdim[5]) ;
  swap_long(&pntr->dime.pixdim[6]) ; swap_long(&pntr->dime.pixdim[7]) ;
  swap_long(&pntr->dime.vox_offset) ; swap_long(&pntr->dime.funused1) ;
  swap_long(&pntr->dime.funused2) ; swap_long(&pntr->dime.cal_max) ;
  swap_long(&pntr->dime.cal_min) ; swap_long(&pntr->dime.compressed) ;
  swap_long(&pntr->dime.verified) ; swap_short(&pntr->dime.dim_un0) ;
  swap_long(&pntr->dime.glmax) ; swap_long(&pntr->dime.glmin) ; }

  swap_long(pntr) unsigned char *pntr;
{ unsigned char b0, b1, b2, b3;

  b0 = *pntr; b1 = *(pntr+1); b2 = *(pntr+2); b3 = *(pntr+3);

  *pntr = b3; *(pntr+1) = b2; *(pntr+2) = b1; *(pntr+3) = b0; }

  swap_short(pntr) unsigned char *pntr;
{ unsigned char b0, b1;

  b0 = *pntr; b1 = *(pntr+1);

  *pntr = b1; *(pntr+1) = b0; }
