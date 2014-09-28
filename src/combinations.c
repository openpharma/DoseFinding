/*
#######################################################################
## This program is Open Source Software: you can redistribute it
## and/or modify it under the terms of the GNU General Public License
## as published by the Free Software Foundation, either version 3 of
## the License, or (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program. If not, see http://www.gnu.org/licenses/.
*/

#include<stdio.h>
#include<stdlib.h>

/* 
   General idea of the algorithm: It is easier to generate all
   possible "border" positions (instead of the actual numbers per
   group). And then count the number of objects between the borders
   (simple differencing).  Example for N=5 and M=3: oo|ooo|o this has
   border positions 2, 5 and number of objects per group is 2, 3, 1
   (2-0, 5-2, 6-5).
 */
void getcomp(int *out, int *work, 
	     int *N, int *B, int *nComp){ 
  int i,j,k,row;
  for(i=0;i<*nComp;i++){
    row = i*(*B+1);
    /* calculate number of obj in each group from borders */
    out[row] = work[0];
    for(j=1;j<*B;j++){
      out[row+j] = work[j]-work[j-1];
    }
    out[row+*B] = *N-work[*B-1];
    /* always increment rightmost number */
    work[*B-1] += 1;
    /* set right numbers to left number */
    for(j= *B-1;j>0;j--){
      if(work[j] == *N+1){
	work[j-1] += 1;
	for(k=j;k<*B;k++){
	  work[k] = work[j-1];
	}
      }
    }
  }
}


