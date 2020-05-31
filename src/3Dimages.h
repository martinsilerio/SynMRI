#ifndef threedimages_h
#define threedimages_h


typedef struct img3d {
   double *intensity;
   int n1;
   int n2;
   int n3;
   int n;
} img3d;

double getVoxelIntensity(img3d*, int, int, int);

int getArrayIndex(img3d*, int, int, int);

void getCoordinates(img3d*, int, int*, int*, int*);

void saveImage(char*, double**, int, int, int);

void saveSlices(img3d*, const char*);

#endif /* threedimages_h*/
