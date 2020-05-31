#include<stdio.h>
#include<math.h>
#include <string.h>
#include "array.h"
#include "3Dimages.h"



double getVoxelIntensity(img3d *img, int i1, int i2, int i3)
{
int i;

i = i3*(img->n1)*(img->n2) +
    i2*(img->n1) +
    i1;

return img->intensity[i];

}

int getArrayIndex(img3d *img, int i1, int i2, int i3)
{
int i;

i = i3*(img->n1)*(img->n2) +
    i2*(img->n1) +
    i1;

return i;

}

void getCoordinates(img3d *img, int index, int *i1, int *i2, int *i3)
{
int res;

*i3 = index/(img->n1)*(img->n2);
res = index%((img->n1)*(img->n2));
*i2 = res/(img->n1);
*i1 = res%(img->n1);

}

void saveImage(char *fileName, double **I, int nr, int nc, int maxim)
   {

   FILE* image;
   int i, j, aux;

   image = fopen(fileName, "w");

   /* here we need a message in case the file can't be opened */
   /*if(!archivo.is_open())
	  {
	  cout << "No se pudo abrir el archivo " << nombreArchivo << endl;
	  getchar();
	  exit(EXIT_FAILURE);
	  } */

   /* first we save the headers, we use the ones GIMP uses */
   fprintf(image, "P2\n");

   fprintf(image, "# CREATOR: GIMP PNM Filter Version 1.1\n");

   /*we save the dimensions*/
   fprintf(image, "%d %d\n", nc, nr);

   /*we save the maximum intensity value*/
   fprintf(image, "%d\n", maxim);

   /*we save the rest of the pixel intensity values*/
   for(i=0;i<nr;i++)
	  for(j=0;j<nc;j++)
		 {
		 aux=round(I[i][j]);
		 if(aux<0)
			aux=0;
		 if(aux>maxim)
			aux=maxim;
         fprintf(image, "%d\n", aux);
		 }

   fclose(image);

   }

void saveSlices(img3d *img, const char* fileName)
{

/* This function saves slides in the middle of each coordinate direction */

int i1_0, i2_0, i3_0;
int i1, i2, i3;
double **I_1, **I_2, **I_3;
double max_1, max_2, max_3;
char *fileN;
char file_1[256] = "s1_";
char file_2[256] = "s2_";
char file_3[256] = "s3_";

MAKE_MATRIX(I_1, img->n2, img->n3);
MAKE_MATRIX(I_2, img->n1, img->n3);
MAKE_MATRIX(I_3, img->n1, img->n2);

i1_0 = img->n1/2;
i2_0 = img->n2/2;
i3_0 = img->n3/2;

max_1 = 0;
for(i2 = 0; i2 < img->n2; i2++)
   for(i3 = 0; i3 < img->n3; i3++)
      {
      I_1[i2][i3] = getVoxelIntensity(img, i1_0, i2, i3);
      if(I_1[i2][i3] > max_1)
         max_1 = round(I_1[i2][i3]);
      }

max_2 = 0;
for(i1 = 0; i1 < img->n1; i1++)
   for(i3 = 0; i3 < img->n3; i3++)
      {
      I_2[i1][i3] = getVoxelIntensity(img, i1, i2_0, i3);
      if(I_2[i1][i3] > max_2)
         max_2 = round(I_2[i1][i3]);
      }

max_3 = 0;
for(i1 = 0; i1 < img->n1; i1++)
   for(i2 = 0; i2 < img->n2; i2++)
      {
      I_3[i1][i2] = getVoxelIntensity(img, i1, i2, i3_0);
      if(I_3[i1][i2] > max_3)
         max_3 = round(I_3[i1][i2]);
      }

fileN = strncat(file_1, fileName, 253);
saveImage(fileN, I_1, img->n2, img->n3,max_1);
fileN = strncat(file_2, fileName, 253);
saveImage(fileN, I_2, img->n1, img->n3,max_2);
fileN = strncat(file_3, fileName, 253);
saveImage(fileN, I_3, img->n1, img->n2,max_3);

FREE_MATRIX(I_1);
FREE_MATRIX(I_2);
FREE_MATRIX(I_3);

}

void saveValues(img3d *img, const char* fileName)
{
FILE *output;
int i;

output = fopen(fileName, "w");
for(i = 0; i < img -> n; i++)
    fprintf(output, "%lf\n", img ->intensity[i]);

fclose(output);

}
