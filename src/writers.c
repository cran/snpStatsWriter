#include <stdio.h>
#include <R.h>

/* Chris Wallace <chris.wallace@cimr.cam.ac.uk> */
/* GPL */

void allelestring(char *A1, char *A2, int idx, char *astr) {
  // simpler version of simple_allelestring.  assumes sep=" ", nullallele="?".
  //  printf("alleles %s, %s, index %s\n",*A1,*A2,*idx);
  char def[] = "? ?";
  strcpy(astr,def);
  switch(idx) {
  case 1 :
    astr[0] = *A1;
    astr[2] = *A1;
    break;
  case 2 :
    astr[0] = *A1;
    astr[2] = *A2;
    break;
  case 3 :
    astr[0] = *A2;
    astr[2] = *A2;
    break;
  default:
    break;
  }
}

void simple_allelestring(char *A1, char *A2, int idx, char *sep, char *nullallele, char *astr) {
  //  printf("alleles %s, %s, index %s\n",*A1,*A2,*idx);
  char def[4];
  def[0]=*nullallele;
  int a1_index=0;
  int a2_index;
  if(strlen(sep)>0) {
    a2_index=2;
    def[1]=*sep;
    def[2]=*nullallele;
    def[3]='\0';
  } else {
    a2_index=1;
    def[1]=*nullallele;
    def[2]='\0';
  }
  strcpy(astr,def);
  switch(idx) {
  case 1 :
    astr[a1_index] = *A1;
    astr[a2_index] = *A1;
    break;
  case 2 :
    astr[a1_index] = *A1;
    astr[a2_index] = *A2;
    break;
  case 3 :
    astr[a1_index] = *A2;
    astr[a2_index] = *A2;
    break;
  default:
    break;
  }
}

void numstring(int idx, char *astr) {
  astr[0] = '4';
  astr[1] = '\0';
  switch(idx) {
  case 1 :
    astr[0] = '1';
    break;
  case 2 :
    astr[0] = '2';
    break;
  case 3 :
    astr[0] = '3';
    break;
  default:
    break;
  }
}

void write_simple(char *x, char **a1, char **a2, int *bp, int *dobp,
		      char **fsep, char **gsep,
		  int *num_coding, // use 1,2,3,4 (4=missing) for alleles, otherwise use characters given in a1, a2
		  char **nullallele, // value for missing alleles
		  int *transpose, // write one snp per line, instead of one sample per line
		      char **file, int *N, int *M,
		  char **rnames, char **cnames, char **eol) {
  int nrow = *N;
  int ncol = *M;
  //  char *eol = "\n";
  //  char *na = "??";
  int i=0, j=0, ij=0;
  FILE * outfile;
  outfile = fopen(*file, "a");
  if (!outfile) {
    error("output file could not be opened");
  }

  char astr[4]; // to hold current allele string

  if(*transpose==1) {
    // format snp a1 a2 a1 a2 etc
    for(j=0; j<ncol; j++) {
      fputs(cnames[j],outfile);
      if(*dobp==1) {
	//	fputs(fsep[0],outfile);
	fprintf(outfile,"%s%d",fsep[0],bp[j]);
      }
      for(i=0, ij=j*nrow; i<nrow; i++, ij++) {
	fputs(fsep[0],outfile);
	if(*num_coding==1) {
	  numstring( (int) x[ij], astr);
	} else {
	  simple_allelestring(a1[j], a2[j], (int) x[ij], gsep[0], nullallele[0], astr);
	}
	fputs(astr,outfile);
      }
      fputs(*eol,outfile);
    }
  } else {
    // each line has format id a1 a2 a1 a2 etc
    for(i=0; i<nrow; i++) {
      fputs(rnames[i], outfile);
      //    fputs(" 1 0 0 M",outfile);
      for(j=0, ij=i; j<ncol; j++, ij+=nrow) {
	fputs(fsep[0],outfile);
	if(*num_coding==1) {
	  numstring( (int) x[ij], astr);
	} else {
	  simple_allelestring(a1[j], a2[j], (int) x[ij], gsep[0], nullallele[0], astr);
	}
	fputs(astr,outfile);
      }
      fputs(*eol,outfile);
    }
  }
  fclose(outfile);
  return;
}

void write_phase(char *x, char **a1, char **a2, char **file, int *N, int *M,
		 char **rnames, char **cnames, 
		 int *print_bp, int *bp, char **eol) { // print a bp line?
  int nrow = *N;
  int ncol = *M;
  char *sep = " ";
  //  char *eol = "\n";
  //  char *na = "??";
  //  char Nch[8];
  FILE * outfile;
  outfile = fopen(*file, "w");
  if (!outfile) {
    error("output file could not be opened");
  }

  // first line - number of samples
  fprintf(outfile, "%d", *N);
  fputs(*eol,outfile);

  // second line - number of snps
  fprintf(outfile, "%d", *M);
  fputs(*eol,outfile);

  // third line - optional - bp
  if( *print_bp == 1) {
    fputs("P",outfile);
    for(int j=0; j<ncol; j++) {
      fprintf(outfile,"%s%d",sep,bp[j]);
    }
    fputs(*eol,outfile);
  }
  
  // fourth line - SSSS - for compatability with phase
  for(int j=0; j<ncol; j++) {
    fputs("S", outfile);
  }
  fputs(*eol,outfile);

  // genotypes - 3 lines per individual
  char astr1[ncol+1], astr2[ncol+1];  // to hold current allele strings
  astr1[ncol] = '\0';
  astr2[ncol] = '\0';
  for(int i=0; i<nrow; i++) {
    fprintf(outfile,"# %s\n",rnames[i]);
    int j, ij;
    for(j=0, ij=i; j<ncol; j++, ij+=nrow) {
      switch(x[ij]) {
      case 0 :
	astr1[j] = '?';
	astr2[j] = '?';
	break;
      case 1 :
	/* if( *num_coding == 1 ) { */
	/*   astr1[j] = '1'; */
	/*   astr2[j] = '1'; */
	/* } else { */
	  astr1[j] = *a1[j];
	  astr2[j] = *a1[j];
	/* } */
	break;
      case 2 :
	/* if( *num_coding == 1 ) { */
	/*   astr1[j] = '1'; */
	/*   astr2[j] = '2'; */
	/* } else { */
	  astr1[j] = *a1[j];
	  astr2[j] = *a2[j];
	/* } */
	break;
      case 3 :
	/* if( *num_coding == 1 ) { */
	/*   astr1[j] = '2'; */
	/*   astr2[j] = '2'; */
	/* } else { */
	  astr1[j] = *a2[j];
	  astr2[j] = *a2[j];
	/* } */
	break;
      }
    }
    fputs(astr1,outfile);
    fputs(*eol,outfile);
    fputs(astr2,outfile);
    fputs(*eol,outfile);
  }
  fclose(outfile);
  return;
}

/* 2.1 Beagle file format */
/*   Beagle input files have a simple format: rows are variables and columns are individuals. */
/* Here is an example of a Beagle file with three individuals and three genotyped markers: */
/* Example 1 - Sample Beagle file */
/* I       id            1001   1001    1002    1002    1003    1003 */
/* A       diabetes      1      1       2       2       2       2 */
/* M       rs2289311     A      G       G       G       A       G */
/* M       rs1248628     T      T       T       C       T       T */
/* M       rs10762764    G      T       T       T       G       T */


void write_beagle(char *x, char **a1, char **a2, int *bp, int *trait,
		  char **gfile, char **mfile, int *N, int *M, int *Ntrait,
		  char **rnames, char **cnames, char **eol) {
  int nrow = *N;
  int ncol = *M;
  int ntrait = *Ntrait;
  char *sep = " ";
  //  char *eol = "\n";
  //  char *na = "??";
  int i=0, j=0, ij=0;
  FILE * outfile_geno;
  FILE * outfile_markers;
  outfile_geno = fopen(*gfile, "w");
  if (!outfile_geno) {
    error("genotype output file could not be opened");
  }
  outfile_markers = fopen(*mfile, "w");
  if (!outfile_markers) {
    error("marker output file could not be opened");
  }

  // markers file
  for(j=0; j<ncol; j++) {
    fprintf(outfile_markers,"%s %d %s %s\n",cnames[j],bp[j],a1[j],a2[j]);
  }

  char astr[4]; // to hold current allele string

  // header
  fprintf(outfile_geno,"I id" );
  for(i=0; i<nrow; i++) {
    fputs(sep,outfile_geno);
    fputs(rnames[i],outfile_geno);
    fputs(sep,outfile_geno);
    fputs(rnames[i],outfile_geno);
  }
  fputs(*eol,outfile_geno);

  // trait?
  if(ntrait>0) {
    fprintf(outfile_geno,"A trait" );
    for(i=0; i<nrow; i++) {
      fprintf(outfile_geno," %d %d",trait[i],trait[i]);
    }
    fputs(*eol,outfile_geno);
  }

// each line has format M rs A1 A2 A1 A2
  for(j=0; j<ncol; j++) {
    fprintf(outfile_geno,"M %s",cnames[j]);
    for (i=0; i<nrow; i++) {
       fputs(sep,outfile_geno);
       allelestring(a1[j], a2[j], (int) x[ij++], astr);
       fputs(astr,outfile_geno);
    }
    fputs(*eol,outfile_geno);
  }

  fclose(outfile_geno);
  fclose(outfile_markers);
  return;
}


void write_mach(char *x, char **a1, char **a2, char **file, 
		char **pedigree, char **member, char **father, char **mother, char **sex,
		int *N, int *M,
		char **rnames, char **cnames, char **eol) {
  int nrow = *N;
  int ncol = *M;
  char *sep = " ";
  //  char *eol = "\n";
  //  char *na = "??";
  int i=0, j=0, ij=0;
  FILE * outfile;
  outfile = fopen(*file, "w");
  if (!outfile) {
    error("output file could not be opened");
    return;
  }

  char astr[4]; // to hold current allele string

// each line has format fam 1 0 0 sex a1 a2 a1 a2 etc
  for(i=0; i<nrow; i++) {
    fputs(pedigree[i], outfile);
    fputs(sep,outfile);
    fputs(member[i], outfile);
    fputs(sep,outfile);
    fputs(father[i], outfile);
    fputs(sep,outfile);
    fputs(mother[i], outfile);
    fputs(sep,outfile);
    fputs(sex[i], outfile);
    fputs(sep,outfile);
    for (j=0, ij=i; j<ncol; j++, ij+=nrow) {
       fputs(sep,outfile);
       allelestring(a1[j], a2[j], (int) x[ij], astr);
       fputs(astr,outfile);
    }
    fputs(*eol,outfile);
  }

  fclose(outfile);
  return;
}

void impute_allelestring(int idx, char *astr) {
  //  printf("alleles %s, %s, index %s\n",*A1,*A2,*idx);
  // write allelestrings as genotype probabilities
  switch(idx) {
  case 1 :
    strcpy(astr,"1 0 0");
    break;
  case 2 :
    strcpy(astr,"0 1 0");
    break;
  case 3 :
    strcpy(astr,"0 0 1");
    break;
  default:
    strcpy(astr,"0 0 0");
    break;
  }
}

void write_impute(char *x, char **a1, char **a2, int *bp, char **file, int *N, int *M,
		  char **rnames, char **cnames, char **snpid, char **eol) {
  int nrow = *N;
  int ncol = *M;
  char *sep = " ";
  //  char *eol = "\n";
  //  char *na = "??";
  int i=0, j=0, ij=0;
  FILE * outfile;
  outfile = fopen(*file, "w");
  if (!outfile) {
    error("output file could not be opened");
    return;
  }

  char astr[6]; // to hold current allele string

// each line has format SNPj rs bp A1 A2 p11 p12 p22 p11 p12 p22 etc
  ij=0;
  for(j=0; j<ncol; j++) {
    fprintf(outfile,"%s %s %d %s %s",snpid[j],cnames[j],bp[j],a1[j],a2[j]);
    for (i=0; i<nrow; i++) {
       fputs(sep,outfile);
       impute_allelestring((int) x[ij++], astr);
       fputs(astr,outfile);
    }
    fputs(*eol,outfile);
  }

  fclose(outfile);
  return;
}

