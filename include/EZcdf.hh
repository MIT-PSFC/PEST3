#ifndef __ezcdfpp__
#define __ezcdfpp__

#include <iostream>
#include <string>
#include "netcdf.h"
#include <stdio.h>

using namespace std;

//namespace EZcdf_use {

  //#include "int2char.h"

//}

const int CDF_MAXDIM = 3; 
const std::string CDF_CHAR = "CHAR";
const std::string CDF_INT = "INT";
const std::string CDF_FLOAT = "FLOAT";
const std::string CDF_DOUBLE = "DOUBLE";



const std::string EZcdf_ErrTag = "--EZcdf Error-- ";
const std::string EZcdf_ErrInFile = EZcdf_ErrTag+" In file ";
const std::string EZcdf_ErrOpn = EZcdf_ErrTag+" Cannot open file ";
const std::string EZcdf_ErrFileNotOpened = EZcdf_ErrTag+" File was not opened ";
const std::string EZcdf_ErrCls = EZcdf_ErrTag+" Cannot close file ";
const std::string EZcdf_ErrType = EZcdf_ErrTag+" Unsupported type ";
const std::string EZcdf_ErrInqType = EZcdf_ErrTag+" Failed to inquire about type for ";
const std::string EZcdf_ErrFindVar = EZcdf_ErrTag+" Cannot find variable ";
const std::string EZcdf_ErrUnknown = EZcdf_ErrTag+" Something strange happened ";
const std::string EZcdf_ErrNotScalar = EZcdf_ErrTag+" Not a scalar ";
const std::string EZcdf_ErrRank = EZcdf_ErrTag+" Rank error for ";
const std::string EZcdf_ErrRankTooLarge = EZcdf_ErrTag+" Rank too large for ";
const std::string EZcdf_ErrSize = EZcdf_ErrTag+" Size error ";
const std::string EZcdf_ErrDefVar = EZcdf_ErrTag+" Undefined variable ";
const std::string EZcdf_ErrGet = EZcdf_ErrTag+" Failed to get ";
const std::string EZcdf_ErrPut = EZcdf_ErrTag+" Failed to put ";
const std::string EZcdf_ErrDefDim = EZcdf_ErrTag+" Failed to define dimension for variable ";
const std::string EZcdf_ErrInq = EZcdf_ErrTag+" Failed to inquire ";
const std::string EZcdf_ErrId = EZcdf_ErrTag+" Failed to get ID for ";
const std::string EZcdf_ErrDimId = EZcdf_ErrTag+" Failed to get dim ID for ";
const std::string EZcdf_ErrDim = EZcdf_ErrTag+" Failed to get dimension length for ";

class EZcdf_Error{
public:
  EZcdf_Error(std::string mes){
    cout << mes.c_str() << endl;
  }
};

/******************************************
  EZcdf: An easy-to-use interface to NetCDF 
  supporting finite sized arrays of 
  depth <= 3
  
  pletzer@pppl.gov September 1999
  ****************************************/

class EZcdf {

private:

  int ncid;
  std::string fname;
  int mode;
  int ier;

  void cdfDefScal(const std::string varnam, std::string xtype);
  nc_type cdfDefType(std::string xtype);


public:

  EZcdf(const std::string fname_in, const std::string m);
  ~EZcdf();

  void cdfDefVar(const std::string varnam, int dimlens[3], std::string xtype);

  void cdfPutVar(const std::string varnam, char *varval);
  void cdfPutVar(const std::string varnam, int *varval);
  void cdfPutVar(const std::string varnam, float *varval);
  void cdfPutVar(const std::string varnam, double *varval);

  void cdfInqVar(const std::string varnam, int dimlens[CDF_MAXDIM], std::string &xtype);

  void cdfGetVar(const std::string varnam, char *varval);
  void cdfGetVar(const std::string varnam, int *varval);
  void cdfGetVar(const std::string varnam, float *varval);
  void cdfGetVar(const std::string varnam, double *varval);


   

};

/* C'tor */

EZcdf::EZcdf(const std::string fname_in, const std::string m="r"){
  fname = fname_in;
  mode = NC_NOCLOBBER;
  if (m=="w" || m=="W") {
    mode = NC_CLOBBER;
    EZcdf::ier = nc_create(fname.c_str(), mode, &ncid);
    if(EZcdf::ier != NC_NOERR) {
      throw EZcdf_Error(EZcdf_ErrOpn+" in write mode (W)");
    }
  } else {
    mode = NC_NOWRITE;
    EZcdf::ier = nc_open(fname.c_str(), mode, &ncid);
    if(EZcdf::ier != NC_NOERR) {
      throw EZcdf_Error(EZcdf_ErrOpn+" in read mode (R)");
    }
  }
}
  
/* D'tor */

EZcdf::~EZcdf(){
  EZcdf::ier = nc_close(ncid);
  if(EZcdf::ier != NC_NOERR) {
    throw EZcdf_Error(EZcdf_ErrCls);
  }
}

/** 
 *  Define Variable. 
 *  Call this prior to writing any variables into the file.
 */

void EZcdf::cdfDefVar(const std::string varnam, int dimlens[CDF_MAXDIM], std::string xtype){

  // check if scalar
  int prod = 1;
  for(int i=0; i<CDF_MAXDIM; i++) prod *= dimlens[i];

  if(prod == 1) {

    cdfDefScal(varnam, xtype);

  } else {

    // if data is array...

    int flag = 0;
    int ndims = 0;
    size_t *dims=NULL; 
    std::string *dnam=NULL;
    int *dimid=NULL;
    int j = 0;

    //EZcdf::ier = nc_redef(ncid);

    // select data type
    nc_type xt = cdfDefType(xtype);

    for(int i=0; i<CDF_MAXDIM; i++){
      if(dimlens[i] > 1 && flag==0){
	ndims = CDF_MAXDIM - i;
	dims = new size_t[ndims];
	dnam = new std::string[ndims];
	dimid = new int[ndims];
	flag = 1;
      }
      if(flag==1){
	char number[2];
	sprintf(number, "%d", CDF_MAXDIM-i);
	dnam[j] =  varnam+"_dim"+std::string(number);
	//dnam[j] =  varnam+"_dim"+std::string(EZcdf_use::Int2Chars(CDF_MAXDIM-i));
	dims[j] = size_t(dimlens[i]);
	EZcdf::ier = nc_def_dim(ncid, dnam[j].c_str(), dims[j], &dimid[j]);
	if(EZcdf::ier != NC_NOERR) {
	  throw EZcdf_Error(EZcdf_ErrDefDim+varnam);
	}
	j++;
      }
    }
  int varid;
  EZcdf::ier = nc_def_var(ncid, varnam.c_str(), xt, ndims, dimid, &varid);
  if(EZcdf::ier != NC_NOERR) throw EZcdf_Error(EZcdf_ErrDefVar+varnam);
  delete[] dims;
  delete[] dnam;
  delete[] dimid;
  }
}


nc_type EZcdf::cdfDefType(std::string xtype){
 // select data type
  nc_type xt;
  std::string xtype_s(xtype);
  if(xtype_s=="CHAR"){
    xt = NC_CHAR;
  }else if(xtype_s=="INT"){
    xt = NC_INT;
  }else if(xtype_s=="R4" || xtype_s=="FLOAT"){
    xt = NC_FLOAT;
  }else if(xtype_s=="R8" || xtype_s=="DOUBLE"){
    xt = NC_DOUBLE;
  }else{
    throw EZcdf_Error(EZcdf_ErrType);
  }
  return xt;
}

void EZcdf::cdfDefScal(const std::string varnam,  std::string xtype){
 // select data type
  nc_type xt = cdfDefType(xtype);
  int dimid;
  int varid;
  EZcdf::ier = nc_def_var(ncid, varnam.c_str(), xt, 0, &dimid, &varid);
}

/**
 *  Put (write) data 
 */

/* CHAR */

void EZcdf::cdfPutVar(const std::string varnam, char *varval){

  // switch to data mode
  EZcdf::ier = nc_enddef(ncid);

  int varid;
  EZcdf::ier = nc_inq_varid (ncid, varnam.c_str(), &varid);
  if(EZcdf::ier != NC_NOERR) {
    throw EZcdf_Error(EZcdf_ErrId+varnam);
  }

  // write data
  EZcdf::ier = nc_put_var_text(ncid, varid, varval);
  if(EZcdf::ier != NC_NOERR) {
    throw EZcdf_Error(EZcdf_ErrPut+varnam);
  }
}

/* INT */

void EZcdf::cdfPutVar(const std::string varnam, int *varval){

  // switch to data mode
  EZcdf::ier = nc_enddef(ncid);

  int varid;
  EZcdf::ier = nc_inq_varid (ncid, varnam.c_str(), &varid);
  if(EZcdf::ier != NC_NOERR) {
    throw EZcdf_Error(EZcdf_ErrId+varnam);
  }

  // write data
  EZcdf::ier = nc_put_var_int(ncid, varid, varval);
  if(EZcdf::ier != NC_NOERR) {
    throw EZcdf_Error(EZcdf_ErrPut+varnam);
  }
}


/* FLOAT */

void EZcdf::cdfPutVar(const std::string varnam, float *varval){

  // switch to data mode
  EZcdf::ier = nc_enddef(ncid);

  int varid;
  EZcdf::ier = nc_inq_varid (ncid, varnam.c_str(), &varid);
  if(EZcdf::ier != NC_NOERR) {
    throw EZcdf_Error(EZcdf_ErrId+varnam);
  }

  // write data
  EZcdf::ier = nc_put_var_float(ncid, varid, varval);
  if(EZcdf::ier != NC_NOERR) {
    throw EZcdf_Error(EZcdf_ErrPut+varnam);
  }
}

/* DOUBLE */

void EZcdf::cdfPutVar(const std::string varnam, double *varval){

  // switch to data mode
  EZcdf::ier = nc_enddef(ncid);

  int varid;
  EZcdf::ier = nc_inq_varid (ncid, varnam.c_str(), &varid);
  if(EZcdf::ier != NC_NOERR) {
    throw EZcdf_Error(EZcdf_ErrId+varnam);
  }

  // write data
  EZcdf::ier = nc_put_var_double(ncid, varid, varval);
  if(EZcdf::ier != NC_NOERR) {
    throw EZcdf_Error(EZcdf_ErrPut+varnam);
  }
}



/**
 *  Inquire about a Variable. 
 *  Return the data shape and type of varnam.
 */

void EZcdf::cdfInqVar(const std::string varnam, int dimlens[CDF_MAXDIM], std::string &xtype){
  int varid;
  EZcdf::ier = nc_inq_varid (ncid, varnam.c_str(), &varid);
  if(EZcdf::ier != NC_NOERR) {
    std::string mes = EZcdf_ErrId+varnam;
    throw EZcdf_Error(mes);
    cout << "in EZcdf::cdfInqVar" << endl;
  } else {

    // get depth
    int ndims;
    EZcdf::ier = nc_inq_varndims(ncid, varid, &ndims);
    if(EZcdf::ier != NC_NOERR) {
      throw EZcdf_Error(EZcdf_ErrRank+varnam);
    }
    if(ndims > CDF_MAXDIM){
      throw EZcdf_Error(EZcdf_ErrRankTooLarge+varnam);
    }

    // get dimensions
    int *ids = new int[ndims];
    EZcdf::ier = nc_inq_vardimid(ncid, varid, ids);
    if(EZcdf::ier != NC_NOERR) {
      throw EZcdf_Error(EZcdf_ErrDimId+varnam);
    }

    for(int i=0; i<CDF_MAXDIM; i++) dimlens[i] = 1; // default is scalar

    for(int i=0; i<ndims; i++){
      size_t length;
      EZcdf::ier = nc_inq_dimlen(ncid, ids[i], &length);
      if(EZcdf::ier != NC_NOERR) {
	throw EZcdf_Error(EZcdf_ErrDim+varnam);
      }
      dimlens[i+CDF_MAXDIM-ndims] = int(length);
    }
    delete[] ids;

    // get data type
    nc_type xt;
    EZcdf::ier = nc_inq_vartype(ncid, varid, &xt);
    if(EZcdf::ier != NC_NOERR) {
      throw EZcdf_Error(EZcdf_ErrInqType+varnam);
  }
    switch (xt){
    case NC_FLOAT:
      xtype = CDF_FLOAT.c_str();
      break;
    case NC_DOUBLE:
      xtype = CDF_DOUBLE.c_str();
      break;
    case NC_INT:
      xtype = CDF_INT.c_str();
      break;
    case NC_CHAR:
      xtype = CDF_CHAR.c_str();
      break;
    default:
      throw  EZcdf_Error(EZcdf_ErrType+varnam);
      break;
    }
  }
}



/**
 * Get (read) data 
 */


/* CHAR */

void EZcdf::cdfGetVar(const std::string varnam, char *varval){
  int varid;

  // switch to data mode
  EZcdf::ier = nc_enddef(ncid);

  // get the variable id
  EZcdf::ier = nc_inq_varid (ncid, varnam.c_str(), &varid);
  if(EZcdf::ier != NC_NOERR) {
    throw EZcdf_Error(EZcdf_ErrId+varnam);
  }
  
  // read the data
  EZcdf::ier = nc_get_var_text(ncid, varid, varval);
  if(EZcdf::ier != NC_NOERR) {
    throw EZcdf_Error(EZcdf_ErrGet+varnam);
  }
}

/* INT */

void EZcdf::cdfGetVar(const std::string varnam, int *varval){
  int varid;

  // switch to data mode
  EZcdf::ier = nc_enddef(ncid);

  // get the variable id
  EZcdf::ier = nc_inq_varid (ncid, varnam.c_str(), &varid);
  if(EZcdf::ier != NC_NOERR) {
    throw EZcdf_Error(EZcdf_ErrId+varnam);
  }
  
  // read the data
  EZcdf::ier = nc_get_var_int(ncid, varid, varval);
  if(EZcdf::ier != NC_NOERR) {
    throw EZcdf_Error(EZcdf_ErrGet+varnam);
  }
}

/* FLOAT */

void EZcdf::cdfGetVar(const std::string varnam, float *varval){
  int varid;

  // switch to data mode
  EZcdf::ier = nc_enddef(ncid);

  // get the variable id
  EZcdf::ier = nc_inq_varid (ncid, varnam.c_str(), &varid);
  if(EZcdf::ier != NC_NOERR) {
    throw EZcdf_Error(EZcdf_ErrId+varnam);
  }
  
  // read the data
  EZcdf::ier = nc_get_var_float(ncid, varid, varval);
  if(EZcdf::ier != NC_NOERR) {
    throw EZcdf_Error(EZcdf_ErrGet+varnam);
  }
}

void EZcdf::cdfGetVar(const std::string varnam, double *varval){
  int varid;

  // switch to data mode
  EZcdf::ier = nc_enddef(ncid);

  // get the variable id
  EZcdf::ier = nc_inq_varid (ncid, varnam.c_str(), &varid);
  if(EZcdf::ier != NC_NOERR) {
    throw EZcdf_Error(EZcdf_ErrId+varnam);
  }
  
  // read the data
  EZcdf::ier = nc_get_var_double(ncid, varid, varval);
  if(EZcdf::ier != NC_NOERR) {
    throw EZcdf_Error(EZcdf_ErrGet+varnam);
  }
}


#endif /* __ezcdfpp__ */


  
