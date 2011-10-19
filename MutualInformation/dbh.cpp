#include "dbh.h"

void swap_long( void  *p) {

  unsigned char *pntr = (unsigned char *)p;
  unsigned char b0, b1, b2, b3;

  b0 = *pntr; 
  b1 = *(pntr+1); 
  b2 = *(pntr+2); 
  b3 = *(pntr+3);

  *pntr = b3; 
  *(pntr+1) = b2; 
  *(pntr+2) = b1; 
  *(pntr+3) = 
    b0; 
};

void swap_short( void *p) { 

  unsigned char *pntr = (unsigned char *)p;
  unsigned char b0, b1;

  b0 = *pntr; 
  b1 = *(pntr+1);

  *pntr = b1; 
  *(pntr+1) = b0; 
};

void swap_hdr( struct dsr *pntr) {
  swap_long(&pntr->hk.sizeof_hdr) ; 
  swap_long(&pntr->hk.extents) ;
  swap_short(&pntr->hk.session_error) ; 
  swap_short(&pntr->dime.dim[0]) ;
  swap_short(&pntr->dime.dim[1]) ; 
  swap_short(&pntr->dime.dim[2]) ;
  swap_short(&pntr->dime.dim[3]) ; 
  swap_short(&pntr->dime.dim[4]) ;
  swap_short(&pntr->dime.dim[5]) ; 
  swap_short(&pntr->dime.dim[6]) ;
  swap_short(&pntr->dime.dim[7]) ; 
  swap_short(&pntr->dime.unused1) ;
  swap_short(&pntr->dime.datatype) ; 
  swap_short(&pntr->dime.bitpix) ;
  swap_long(&pntr->dime.pixdim[0]) ;
  swap_long(&pntr->dime.pixdim[1]) ;
  swap_long(&pntr->dime.pixdim[2]) ; 
  swap_long(&pntr->dime.pixdim[3]) ;
  swap_long(&pntr->dime.pixdim[4]) ; 
  swap_long(&pntr->dime.pixdim[5]) ;
  swap_long(&pntr->dime.pixdim[6]) ; 
  swap_long(&pntr->dime.pixdim[7]) ;
  swap_long(&pntr->dime.vox_offset) ; 
  swap_long(&pntr->dime.funused1) ;
  swap_long(&pntr->dime.funused2) ; 
  swap_long(&pntr->dime.cal_max) ;
  swap_long(&pntr->dime.cal_min) ; 
  swap_long(&pntr->dime.compressed) ;
  swap_long(&pntr->dime.verified) ;
  swap_short(&pntr->dime.dim_un0) ;
  swap_long(&pntr->dime.glmax) ; 
  swap_long(&pntr->dime.glmin) ; 
};


