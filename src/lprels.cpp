/*============================================================================
    Copyright 2006 William Hart    

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

===============================================================================*/

/* 
   This code has been adapted for FLINT from mpqs.c in the Pari/GP package. 
   See http://pari.math.u-bordeaux.fr/
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <gmp.h>
#include <unistd.h>
#include "lprels.h"

#define min_bufspace 120UL  /* use new buffer when < min_bufspace left */
#define buflist_size 4096UL /* size of list-of-buffers blocks */
#define sort_table_size 100000UL

/*********************************************************************

    File based large prime routines
    
*********************************************************************/

/* 
    Concatenates a filename and directory name to give a full path 
*/

char * get_filename(const char *dir, const char *s)
{
  char *buf = (char *) malloc(strlen(dir) + strlen(s) + 2);
#if defined(__EMX__) || defined(WINCE)
  sprintf(buf, "%s\\%s", dir,s);
#else
  sprintf(buf, "%s/%s", dir,s);
#endif
  return buf;
}

char * unique_filename(const char *s)
{
  char *buf, suf[64];
  size_t lsuf;

  sprintf(suf,".%ld.%ld", (long)getuid(), (long)getpid());

  lsuf = strlen(suf);
  /* room for s + suffix '\0' */
  buf = (char*) malloc(8 + lsuf + 1);
  
  sprintf(buf, "%.8s%s", s, suf);
  return buf;
}


FILE * flint_fopen(const char * name, const char * mode)
{
#if defined(WINCE) || defined(macintosh)
  const char * tmp_dir = NULL;
#else
  const char * tmp_dir = getenv("TMPDIR");
#endif
  if (tmp_dir == NULL) tmp_dir = "./";
  FILE * temp_file = fopen(get_filename(tmp_dir,unique_filename(name)),mode);
  if (!temp_file)
  {
     printf("Unable to open temporary file\n");
     abort();
  }
  return temp_file;
}

/* 
    Compares two large prime relations according to their first 
    element (the large prime). Used by qsort.
*/

int relations_cmp(const void *a, const void *b)
{
  char **sa = (char**) a;
  char **sb = (char**) b;
  long qa = strtol(*sa, NULL, 10);
  long qb = strtol(*sb, NULL, 10);
  
  if (qa < qb) return -1;
  else if (qa > qb) return 1;
  else return strcmp(*sa, *sb);
}


/*
    Writes the given string to the given file and aborts upon error
*/

void flint_fputs(const char *s, FILE *file)
{
  if (fputs(s, file) < 0)
  {
      printf("Error whilst writing to large prime file!");
      abort();
  }
}

/* 
   Given a file "filename" containing full or large prime relations,
   rearrange the file so that relations are sorted by their first elements.
   Works in memory, discards duplicate lines, and overwrites the original 
   file. Returns the number of relations after sorting and discarding.
*/

long sort_lp_file(const char *filename)
{
  FILE *TMP;
  char *old_s, *buf, *cur_line;
  char **s_table, **sort_table, **buflist, **buflist_head;
  long i, j, count;
  size_t length, bufspace;
  
  buflist_head = (char**) malloc(buflist_size * sizeof(char*));
  buflist = buflist_head;
  
  *buflist++ = NULL; /* flag this as last and only buflist block */
  
  TMP = flint_fopen(filename, "r");
  
  /* allocate first buffer and read first line, if any, into it */
  buf = (char*) malloc(MPQS_STRING_LENGTH * sizeof(char));
  cur_line = buf;
  bufspace = MPQS_STRING_LENGTH;

  if (fgets(cur_line, bufspace, TMP) == NULL)
  { /* file empty */
    free(buf); free(buflist_head); 
    fclose(TMP);
    return 0;
  }
  /* enter first buffer into buflist */
  *buflist++ = buf; /* can't overflow the buflist block */
  length = strlen(cur_line) + 1; /* count the \0 byte as well */
  bufspace -= length;

  s_table = (char**) malloc(sort_table_size*sizeof(char*));
  sort_table = s_table+sort_table_size;
  
  /* at start of loop, one line from the file is sitting in cur_line inside buf,
     the next will go into cur_line + length, and there's room for bufspace
     further characters in buf. The loop reads another line if one exists, and
     if this overruns the current buffer, it allocates a fresh one --GN */
   
  for (i = 0, sort_table--; /* until end of file */; i++, sort_table--)
  { 
       
    *sort_table = cur_line;
    cur_line += length;

    /* if little room is left, allocate a fresh buffer before attempting to
     * read a line, and remember to free it if no further line is forthcoming.
     * This avoids some copying of partial lines --GN */
    if (bufspace < min_bufspace)
    {
      buf = (char*) malloc(MPQS_STRING_LENGTH * sizeof(char));
      cur_line = buf;
      bufspace = MPQS_STRING_LENGTH;
      if (fgets(cur_line, bufspace, TMP) == NULL) { free(buf); break; }

      if (buflist - buflist_head >= buflist_size) abort();
      
      /* remember buffer for later deallocation */
      *buflist++ = buf;
      length = strlen(cur_line) + 1;
      bufspace -= length; continue;
    }

    /* normal case:  try fitting another line into the current buffer */
    if (fgets(cur_line, bufspace, TMP) == NULL) break; /* none exists */
    length = strlen(cur_line) + 1;
    bufspace -= length;

    /* check whether we got the entire line or only part of it */
    if (bufspace == 0 && cur_line[length-2] != '\n')
    {
      size_t lg1;
      buf = (char*) malloc(MPQS_STRING_LENGTH * sizeof(char));
      if (buflist - buflist_head >= buflist_size) abort();
      /* remember buffer for later deallocation */
      *buflist++ = buf;

      /* copy what we've got to the new buffer */
      (void)strcpy(buf, cur_line); /* cannot overflow */
      cur_line = buf + length - 1; /* point at the \0 byte */
      bufspace = MPQS_STRING_LENGTH - length + 1;
      
      /* read remainder of line */
      if (fgets(cur_line, bufspace, TMP) == NULL)
      {
         printf("MPQS: relations file truncated?!\n");
         abort();
      }
      lg1 = strlen(cur_line);
      length += lg1; /* we already counted the \0 once */
      bufspace -= (lg1 + 1); /* but here we must take it into account */
      cur_line = buf; /* back up to the beginning of the line */
    }
  } /* for */

  fclose(TMP);

  /* sort the whole lot in place by swapping pointers */
  qsort(sort_table, i, sizeof(char*), relations_cmp);

  /* copy results back to the original file, skipping exact duplicates */
  TMP = flint_fopen(filename, "w");
  old_s = sort_table[0];
  flint_fputs(sort_table[0], TMP);
  count = 1;
  for(j = 1; j < i; j++)
  {
    if (strcmp(old_s, sort_table[j]))
    {
      flint_fputs(sort_table[j], TMP);
      count++;
    }
    old_s = sort_table[j];
  }
  fflush(TMP);
  fclose(TMP);
  /* deallocate buffers */  
  while (*--buflist)
  {
     if (buflist != buflist_head)
        free(*buflist);   /* free a buffer */
  }
  free(buflist_head); 
  free(s_table);
  return count;
}

/* 
   Appends contents of file fp1 to fp (auxiliary routine for merge sort) and
   returns number of lines copied. Closes fp afterwards.
*/

long append_file(FILE *fp, FILE *fp1)
{
  char line[MPQS_STRING_LENGTH];
  long c = 0;
  while (fgets(line, MPQS_STRING_LENGTH, fp1)) { flint_fputs(line, fp); c++; }
  if (fflush(fp)) 
  {
     printf("Error while flushing file.\n");
     abort();                
  }                
  fclose(fp); return c;
}

/* 
   Merge-sort on the files LPREL and LPNEW; assumes that LPREL and LPNEW are
   already sorted. Creates/truncates the TMP file, writes result to it and
   closes it (via append_file()). Instead of LPREL, LPNEW we may also call
   this with FREL, FNEW. In the latter case COMB should be NULL (and we
   return the count of all full relations), in the former case it should be
   non-NULL (and we return the count of frels we expect to be able to combine 
   out of the present lprels). If COMB is non-NULL, the combinable lprels 
   are written out to this separate file.
   
   We retain only one occurrence of each large prime in TMP (i.e. in the
   future LPREL file). --GN 
*/

#define swap_lines() { char *line_tmp;\
  line_tmp = line_new_old; \
  line_new_old = line_new; \
  line_new = line_tmp; }

long mergesort_lp_file_internal(FILE *LPREL, FILE *LPNEW, FILE *COMB, FILE *TMP)
{
  char line1[MPQS_STRING_LENGTH], line2[MPQS_STRING_LENGTH];
  char line[MPQS_STRING_LENGTH];
  char *line_new = line1, *line_new_old = line2;
  long q_new, q_new_old = -1, q, i = 0, c = 0;
  long comb_in_progress;

  if ( !fgets(line_new, MPQS_STRING_LENGTH, LPNEW) )
  { 
    /* LPNEW is empty: copy LPREL to TMP. Could be done by a rename if we
       didn't want to count the lines (again)... however, this case will not
       normally happen */
    i = append_file(TMP, LPREL);
    return COMB ? 0 : i;
  }
  /* we now have a line_new from LPNEW */

  if (!fgets(line, MPQS_STRING_LENGTH, LPREL))
  { 
    /* LPREL is empty: copy LPNEW to TMP... almost. */
    flint_fputs(line_new, TMP);
    
    if (!COMB)
    { 
      /* full relations mode */
      i = append_file(TMP, LPNEW);
      return i + 1;
    }

    /* LP mode:  check for combinable relations */
    q_new_old = atol(line_new);
    /* we need to retain a copy of the old line just for a moment, because we
       may yet have to write it to COMB. Do this by swapping the two buffers */
    swap_lines();
    comb_in_progress = 0;
    i = 0;

    while (fgets(line_new, MPQS_STRING_LENGTH, LPNEW))
    {
      q_new = atol(line_new);
      if (q_new_old == q_new)
      { 
        /* found combinables, check whether we're already busy on this
           particular large prime */
        if (!comb_in_progress)
        { 
          /* if not, write first line to COMB, creating and opening the
             file first if it isn't open yet */
          flint_fputs(line_new_old, COMB);
          comb_in_progress = 1;
        }
        /* in any case, write the current line, and count it */
        flint_fputs(line_new, COMB);
        i++;
      }
      else
      { 
        /* not combinable */
        q_new_old = q_new;
        comb_in_progress = 0;
        /* and dump it to the TMP file */
        flint_fputs(line_new, TMP);
        /* and stash it away for a moment */
        swap_lines();
        comb_in_progress = 0;
      }
    } /* while */
    fclose(TMP); return i;
  }

  /* normal case: both LPNEW and LPREL are not empty */
  q_new = atol(line_new);
  q = atol(line);

  for(;;)
  { 
    /* main merging loop */
    i = comb_in_progress = 0;

    /* first the harder case:  let LPNEW catch up with LPREL, and possibly
       overtake it, checking for combinables coming from LPNEW alone */
    while (q > q_new)
    {
      if (!COMB || !comb_in_progress) flint_fputs(line_new, TMP);
      if (!COMB) c++; /* in FREL mode, count lines written */
      else if (!comb_in_progress)
      {
        q_new_old = q_new;
        swap_lines();
      }
      if (!fgets(line_new, MPQS_STRING_LENGTH, LPNEW))
      {
        flint_fputs(line, TMP);
        if (!COMB) c++; else c += i;
        i = append_file(TMP, LPREL);
        return COMB ? c : c + i;
      }
      q_new = atol(line_new);
      if (!COMB) continue;

      /* LP mode only: */
      if (q_new_old != q_new) /* not combinable */
        comb_in_progress = 0; /* next loop will deal with it, or loop may end */
      else
      { 
        /* found combinables, check whether we're already busy on this
           large prime */
        if (!comb_in_progress)
        {
          flint_fputs(line_new_old, COMB);
          comb_in_progress = 1;
        }
        /* in any case, write the current line, and count it */
        flint_fputs(line_new, COMB);
        i++;
      }
    } /* while q > q_new */

    /* q <= q_new */

    if (COMB) c += i;    /* accumulate count of combinables */
    i = 0;               /* and clear it */
    comb_in_progress = 0;/* redundant */

    /* now let LPREL catch up with LPNEW, and possibly overtake it */
    while (q < q_new)
    {
      flint_fputs(line, TMP);
      if (!COMB) c++;
      if (!fgets(line, MPQS_STRING_LENGTH, LPREL))
      {
        flint_fputs(line_new, TMP);
        i = append_file(TMP, LPNEW);
        return COMB ? c : c + i + 1;
      }
      else
        q = atol(line);
    }

    /* q >= q_new */

    /* Finally, it may happen that q == q_new, indicating combinables whose
       large prime is already in LPREL, and appears now one or more times in
       LPNEW. Thus in this sub-loop we advance LPNEW. The `line' from LPREL is
       left alone, and will be written to TMP the next time around the main for
       loop; we only write it to COMB here -- unless all we find is an exact
       duplicate of the line we already have, that is. (There can be at most
       one such, and if so it is simply discarded.) */
    while (q == q_new)
    {
      if (!strcmp(line_new, line))
      { 
        /* duplicate -- move right ahead to the next LPNEW line */
        ;/* do nothing here */
      }
      else if (!COMB)
      { /* full relations mode: write line_new out first, keep line */
        flint_fputs(line_new, TMP);
        c++;
      }
      else
      { 
        /* LP mode, and combinable relation */
        if (!comb_in_progress)
        {
          flint_fputs(line, COMB);
          comb_in_progress = 1;
        }
        flint_fputs(line_new, COMB);
        i++;
      }
      /* NB comb_in_progress is cleared by q_new becoming bigger than q, thus
         the current while loop terminating, the next time through the main for
         loop */

      /* common ending: get another line_new, if any */
      if (!fgets(line_new, MPQS_STRING_LENGTH, LPNEW))
      {
        flint_fputs(line, TMP);
        if (!COMB) c++; else c += i;
        i = append_file(TMP, LPREL);
        return COMB ? c : c + i;
      }
      else
        q_new = atol(line_new);
    } /* while */

    if (COMB) c += i; /* accumulate count of combinables */
  }
}

/* 
   Perform mergesort of large prime files
*/

long mergesort_lp_file(const char *REL_str, const char *NEW_str, const char *TMP_str, FILE *COMB)
{
  FILE *NEW = flint_fopen(NEW_str, "r");
  
#if defined(WINCE) || defined(macintosh)
  const char * tmp_dir = NULL;
#else
  const char * tmp_dir = getenv("TMPDIR");
#endif
  if (tmp_dir == NULL) tmp_dir = "./";
  char * TMP_name = get_filename(tmp_dir,unique_filename(TMP_str));
  char * REL_name = get_filename(tmp_dir,unique_filename(REL_str));
  FILE * TMP = fopen(TMP_name,"w");
  FILE * REL = fopen(REL_name,"r");
  if ((!TMP) || (!REL))
  {
     printf("Unable to open temporary file\n");
     abort();
  }
  
  long tp = mergesort_lp_file_internal(REL, NEW, COMB, TMP);
  fclose(REL);
  fclose(NEW);
  
  if (rename(TMP_name,REL_name))
  {
     printf("Cannot rename file %s to %s", TMP_str, REL_str);
     abort();
  } 
  return tp;
}

void read_matrix(unsigned long ** relations, FILE *FREL, la_col_t* colarray, unsigned long * relsFound, unsigned long relSought, mpz_t * XArr, mpz_t n, unsigned long * factorBase)
{
  long e, p;
  char buf[MPQS_STRING_LENGTH], *s;
  //char buf2[MPQS_STRING_LENGTH];
  unsigned long numfactors;
  mpz_t test1, test2;
  mpz_init(test1);
  mpz_init(test2);
   

  if (ftell(FREL) < 0)
  {
     printf("Error on full relations file\n");
     abort();
  }
  while ((fgets(buf, MPQS_STRING_LENGTH, FREL)) && ((*relsFound) < relSought))
  {
    numfactors = 0;
    gmp_sscanf(buf,"%Zd",XArr[*relsFound]);
    s = strchr(buf, ':') + 2;
    s = strtok(s, " \n");
    while (s != NULL)
    {
      e = atol(s); if (!e) break;
      s = strtok(NULL, " \n");
      p = atol(s);
      if (e & 1) xorColEntry(colarray,*relsFound,p);
      for (long i = 0; i < e; i++) relations[*relsFound][++numfactors] = p;
      s = strtok(NULL, " \n");
    }
    relations[*relsFound][0] = numfactors;
    
    mpz_set_ui(test1,1);
    for (unsigned long i = 1; i<=relations[*relsFound][0]; i++)
    {
       mpz_mul_ui(test1,test1,factorBase[relations[*relsFound][i]]);
       if (i%30 == 0) mpz_mod(test1,test1,n);
    }
    mpz_mod(test1,test1,n);
    mpz_mul(test2,XArr[*relsFound],XArr[*relsFound]);
    mpz_mod(test2,test2,n);
    if (mpz_cmp(test1,test2)!=0)
    {
       mpz_add(test1,test1,test2);
       if (mpz_cmp(test1,n)!=0) 
       {
          clearCol(colarray,*relsFound);
       }
       else (*relsFound)++;
    } else (*relsFound)++;  
  }
  
  mpz_clear(test1);
  mpz_clear(test2);
   
  return;
}

/*********************************************************************

    Routines for writing relations as strings
    
*********************************************************************/

/* 
    Writes a factor pi^ei into a string as " ei pi" 
*/

void add_factor(char **last, unsigned long ei, unsigned long pi) 
{
  sprintf(*last, " %ld %ld", ei, pi);
  *last += strlen(*last);
}

/*
   Concatenate " 0" to string 
*/

void add_0(char **last) 
{
  char *s = *last;
  *s++ = ' ';
  *s++ = '0';
  *s++ = 0; *last = s;
}

/*********************************************************************

    Large prime relation combining
    
*********************************************************************/

/*
   Add to an array of unsigned longs the exponents from a relation string
*/

void set_exponents(unsigned long *ei, char *r)
{
  char *s, b[MPQS_STRING_LENGTH];
  long e;

  strcpy(b, r);
  s = strtok(b, " \n");
  while (s != NULL)
  {
    e = atol(s); if (!e) break;
    s = strtok(NULL, " \n");
    ei[atol(s)] += e;
    s = strtok(NULL, " \n");
  }
}

/* Writes an lp_entry from a string */

void set_lp_entry(mpqs_lp_entry *e, char *buf)
{
  char *s1, *s2;
  s1 = buf; s2 = strchr(s1, ' '); *s2 = '\0';
  e->q = atol(s1);
  s1 = s2 + 3; s2 = strchr(s1, ' '); *s2 = '\0';
  strcpy(e->Y, s1);
  s1 = s2 + 3; s2 = strchr(s1, '\n'); *s2 = '\0';
  strcpy(e->E, s1);
}

/* 
   Combines the large prime relations in COMB to full relations in FNEW.
   FNEW is assumed to be open for writing / appending. 
*/

long combine_large_primes(unsigned long numPrimes,
                          FILE *COMB, FILE *FNEW, mpz_t N, mpz_t factor)
{
  char new_relation[MPQS_STRING_LENGTH], buf[MPQS_STRING_LENGTH];
  mpqs_lp_entry e[2]; /* we'll use the two alternatingly */
  unsigned long *ei; 
  long ei_size = numPrimes;
  long old_q;
  mpz_t inv_q, Y1, Y2, new_Y, new_Y1;
  mpz_init(inv_q);mpz_init(Y1);mpz_init(Y2);mpz_init(new_Y);mpz_init(new_Y1);
  long i, l, c = 0;

  if (!fgets(buf, MPQS_STRING_LENGTH, COMB)) return 0; /* should not happen */

  ei = (unsigned long *) malloc(sizeof(unsigned long)*ei_size);
  
  /* put first lp relation in row 0 of e */
  set_lp_entry(&e[0], buf);

  i = 1; /* second relation will go into row 1 */
  old_q = e[0].q;
  mpz_set_ui(inv_q, old_q);
  while (!mpz_invert(inv_q, inv_q, N)) /* can happen */
  {
    /* We have found a factor. It could be N when N is quite small;  
       or we might just have found a divisor by sheer luck. */
    mpz_gcd_ui(inv_q, N, old_q);
    if (!mpz_cmp(inv_q, N)) /* pity */
    {
      if (!fgets(buf, MPQS_STRING_LENGTH, COMB)) { return 0; }
      set_lp_entry(&e[0], buf);
      old_q = e[0].q; 
      mpz_set_ui(inv_q, old_q);
      continue;
    }
    mpz_set(factor, inv_q);
    free(ei);
    return c;
  }
  gmp_sscanf(e[0].Y, "%Zd", Y1);
  
  while (fgets(buf, MPQS_STRING_LENGTH, COMB))
  {
    set_lp_entry(&e[i], buf);
    if (e[i].q != old_q)
    {
      /* switch to combining a new bunch, swapping the rows */
      old_q = e[i].q;
      mpz_set_ui(inv_q, old_q);
      while (!mpz_invert(inv_q, inv_q, N)) /* can happen */
      {
        mpz_gcd_ui(inv_q, N, old_q);
        if (!mpz_cmp(inv_q, N)) /* pity */
        {
          old_q = -1; /* sentinel */
          continue; /* discard this combination */
        }
        mpz_set(factor, inv_q);
        free(ei);
        return c;
      }
      gmp_sscanf(e[i].Y, "%Zd", Y1);
      i = 1 - i; /* subsequent relations go to other row */
      continue;
    }
    /* count and combine the two we've got, and continue in the same row */
    memset((void *)ei, 0, ei_size * sizeof(long));
    set_exponents(ei, e[0].E);
    set_exponents(ei, e[1].E);
    gmp_sscanf(e[i].Y, "%Zd", Y2);
    
    if (mpz_cmpabs(Y1,Y2)!=0)
    {
       c++;
       mpz_mul(new_Y, Y1, Y2);
       mpz_mul(new_Y, new_Y, inv_q);
       mpz_mod(new_Y, new_Y, N);
    
       mpz_sub(new_Y1, N, new_Y);
       if (mpz_cmpabs(new_Y1, new_Y) < 0) mpz_set(new_Y, new_Y1);
    
       gmp_sprintf(buf, "%Zd\0", new_Y); 
       strcpy(new_relation, buf);
       strcat(new_relation, " :");
       for (l = 0; l < ei_size; l++)
          if (ei[l])
          {
             sprintf(buf, " %ld %ld", ei[l], l);
             strcat(new_relation, buf);
          }
       strcat(new_relation, " 0");
       strcat(new_relation, "\n");

       flint_fputs(new_relation, FNEW);
    }
  } /* while */

  free(ei);
  mpz_clear(inv_q);mpz_clear(Y1);mpz_clear(Y2);mpz_clear(new_Y);mpz_clear(new_Y1);

  return c;
}


