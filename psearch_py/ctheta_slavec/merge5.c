/*
 *     File: merge5.c
 * Language: C 
 *   Author: Kenneth J. Mighell
 *  Version: 0.3.2  2018MAY06
 */

#include<stdio.h>

// Sort data vector in ASCENDING order

#define LENI 10
#define DIMI (LENI)

typedef double KEY_T;
typedef long   POS_T;
typedef long   INDEX_T;

/*
 * Merge functions merges the two sorted parts. 
 * Sorted parts will be from [left, mid] and [mid+1, right].
 */
void Merge
(
 KEY_T dataDA[], 
 KEY_T tmpDA[], 
 POS_T left, 
 POS_T mid, 
 POS_T right
 )
{
  /*
   * Use tmpDA to store the new sorted part
   */
  POS_T pos=0;
  POS_T lpos = left;
  POS_T rpos = mid + 1;
  while(lpos <= mid && rpos <= right)
    {
      if(dataDA[lpos] < dataDA[rpos])
	{
	  tmpDA[pos++] = dataDA[lpos++];
	}
      else
	{
	  tmpDA[pos++] = dataDA[rpos++];
	}
    }
  while(lpos <= mid)  
    tmpDA[pos++] = dataDA[lpos++];
  while(rpos <= right)
    tmpDA[pos++] = dataDA[rpos++];
  int j;
  /*
   * Copy back the sorted tmpDA to the original dataDA 
   */
  for(j = 0;j < pos; j++)
    {
      dataDA[j+left] = tmpDA[j];
    }
  return;
}

/*
 * Merge functions merges the two sorted parts. 
 * Sorted parts will be from [left, mid] and [mid+1, right].
 */
void MergeIndex
(
 KEY_T dataDA[], 
 INDEX_T indexLA[],
 INDEX_T tmpLA[],
 POS_T left, 
 POS_T mid, 
 POS_T right
 )
{
  /*
   * Use tmpLA to store the new sorted part
   */
  POS_T pos = 0;
  POS_T lpos = left;
  POS_T rpos = mid + 1;
  int j;
  while(lpos <= mid && rpos <= right)
    if(dataDA[indexLA[lpos]] < dataDA[indexLA[rpos]])
      tmpLA[pos++] = indexLA[lpos++];
    else 
      tmpLA[pos++] = indexLA[rpos++];
  while(lpos <= mid)  
    tmpLA[pos++] = indexLA[lpos++];
  while(rpos <= right)
    tmpLA[pos++] = indexLA[rpos++];
  /*
   * Copy back the sorted tmpLA to the original indexLA 
   */
  for(j = 0;j < pos; j++)
    indexLA[j+left] = tmpLA[j];
  return;
}

/*
 * This function sorts the dataDA in the range [left,right].
 * That is from index left to index right inclusive
 */
void MergeSort
(
 KEY_T dataDA[], 
 KEY_T tmpDA[], 
 POS_T left, 
 POS_T right
 )
{
  POS_T mid = (left+right)/2;
  /* 
   * We have to sort only when left<right because 
   * when left=right it is already sorted
   */
  if(left<right)
    {
      /* Sort the left part */
      MergeSort(dataDA,tmpDA,left,mid);
      /* Sort the right part */
      MergeSort(dataDA,tmpDA,mid+1,right);
      /* Merge the two sorted parts */
      Merge(dataDA,tmpDA,left,mid,right);
    }
}

/*
 * This function sorts the dataDA in the range [left,right].
 * That is from index left to index right inclusive
 */
void MergeSortIndex
(
 KEY_T dataDA[], 
 INDEX_T indexLA[], 
 INDEX_T tmpLA[], 
 POS_T left, 
 POS_T right
 )
{
  POS_T mid = (left+right)/2;
  /* 
   * We have to sort only when left<right because 
   * when left=right it is already sorted
   */
  if(left<right)
    {
      /* Sort the left part */
      MergeSortIndex(dataDA,indexLA,tmpLA,left,mid);
      /* Sort the right part */
      MergeSortIndex(dataDA,indexLA,tmpLA,mid+1,right);
      /* Merge the two sorted parts */
      MergeIndex(dataDA,indexLA,tmpLA,left,mid,right);
    }
}

#undef INDEXED
#define INDEXED

#define TESTING
#undef TESTING
#ifdef TESTING
int main()
{
#ifdef INDEXED
#warning This is the INDEXED version
#endif
#ifndef INDEXED
#warning This is the NONINDEXED version
#endif
  KEY_T dataDA[DIMI];
#ifdef INDEXED
  INDEX_T indexLA[DIMI];
  INDEX_T tmpLA[DIMI];
#else
  KEY_T tmpDA[DIMI];
#endif
  int number_of_elements = DIMI;
  int i;
  int j;
  int nI;
  nI = LENI;
#ifndef INDEXED
  for( i=0; i<nI; i++ )
    {
      double x;
      x = (nI-1) - i;
      x += (i+1)*0.0001;
      dataDA[i] = x;
    }
  printf("\nRaw data: j,dataDA[j]\n");
  for(j = 0;j < number_of_elements;j++)
    {
      printf("%d %10.4f\n",j,dataDA[j]);
    }
  /* Calling this functions sorts the dataDA */
  MergeSort(dataDA,tmpDA,0,number_of_elements-1); 
  printf("\nSorted data:  j,dataDA[j]\n");
  for(j = 0;j < number_of_elements;j++)
    {
      printf("%d %10.4f\n",j,dataDA[j]);
    }
#endif
#ifdef INDEXED
  for( i=0; i<nI; i++ )
    {
      double x;
      x = (nI-1) - i;
      x += (i+1)*0.0001;
      dataDA[i] = x;
      indexLA[i] = i;
    }
  printf("\nRaw data:  j,dataDA[j],indexLA[j]\n");
  for(j = 0;j < number_of_elements;j++)
    {
      printf("%d %10.4f (%ld)\n",j,dataDA[j],indexLA[j]);
    }
  MergeSortIndex(dataDA,indexLA,tmpLA,0,number_of_elements-1); 
  printf("\nSorted data:  j,dataDA[j],dataDA[indexLA[j]],indexLA[j]\n");
  for(j = 0;j < number_of_elements;j++)
    {
      printf("%d %10.4f --> %10.4f (%ld) \n",j,dataDA[j],dataDA[indexLA[j]],indexLA[j]);
    }
#endif
  printf("\n\nBye!\n");
  return 0;
}
#endif 
#undef LENI
#undef DIMI
//EOF
