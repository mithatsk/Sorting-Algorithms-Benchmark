#include<stdio.h>
#include <time.h>
#include <stdlib.h> 
#include <sys/time.h>
#include<math.h>
void swap(int *,int *);
void quick_sort(int *,int);
int * random_array(int,int,int);
static double clockms (void);
double FindMax(double *,int);
double FindMin(double *,int);
double FindAverage(double *,int);

void quick_sort (int *a, int n) {
    int i, j, p, t;
    if (n < 2)
        return;
    p = a[n / 2];
    for (i = 0, j = n - 1;; i++, j--) {
        while (a[i] < p)
            i++;
        while (p < a[j])
            j--;
        if (i >= j)
            break;
        t = a[i];
        a[i] = a[j];
        a[j] = t;
    }
    quick_sort(a, i);
    quick_sort(a + i, n - i);
}

static void reverse_sort(int number[],int size){
	int i,j,a;
	for (i = 0; i < size; ++i)
    {
        for (j = i + 1; j < size; ++j)
        {
            if (number[i] < number[j])
            {
                a = number[i];
                number[i] = number[j];
                number[j] = a;
            }
        }
	
}}

static void merge (int * arr, int * tmp, size_t beg, size_t mid, size_t end)
{
  size_t ii = beg;
  size_t jj = mid;
  size_t kk = beg;
  
  while (ii < mid && jj < end) {
    if (arr[ii] < arr[jj]) {
      tmp[kk++] = arr[ii++];
    }
    else {
      tmp[kk++] = arr[jj++];
    }
  }
  
  while (ii < mid) {
    tmp[kk++] = arr[ii++];
  }
  while (jj < end) {
    tmp[kk++] = arr[jj++];
  }
  
  for (ii = beg; ii < end; ++ii) {
    arr[ii] = tmp[ii];
  }
}


static void msort_rec (int * arr, int * tmp, int beg, int end)
{
  size_t len, mid;
  
  len = end - beg;
  if (1 >= len) {
    return;
  }
  
  mid = beg + len / 2;
  msort_rec (arr, tmp, beg, mid);
  msort_rec (arr, tmp, mid, end);
  merge (arr, tmp, beg, mid, end);
}


static void merge_sort (int * arr, int len)
{
  int * tmp =(int*) malloc (len * sizeof(int));
  if (NULL == tmp) {
   return;
  }
  msort_rec (arr, tmp, 0, len);
  free (tmp);
}


void swap(int *x, int *y)
{
   int temp;
 
   temp = *y;
   *y   = *x;
   *x   = temp;   
}

int * random_array(int maxval,int minval,int len){
    if(maxval<minval){
        swap(&maxval,&minval);
    }
    srand (time(NULL));
    int * rArray;
    int i;
    
    rArray = (int *)malloc(sizeof(int)*len);
    if(rArray==NULL)
        return 0;
    for(i=0;i<len;i++){
        rArray[i] = rand() % (maxval + 1 - minval) + minval;
    }
    return rArray;
}

double clockms (void)
{
  static struct timeval t0 = { 0, 0 };
  struct timeval t1;
  
  if (0 == t0.tv_sec) {
    if (0 != gettimeofday (&t0, NULL)) {
      return 0;
    }
  }
  if (0 != gettimeofday (&t1, NULL)) {
    return 0;
  }
  
  return (1e3 * (t1.tv_sec - t0.tv_sec) + 1e-3 * (t1.tv_usec - t0.tv_usec))/1000;
}

double FindMax(double * arr,int len){
	int aa;
	double max=-1;
	for(aa=0;aa<len;aa++){
		if(max<arr[aa]){
			max=arr[aa];
		}
	}
	return max;
}

double FindMin(double * arr,int len){
	int bb;
	double min=9999;
	for(bb=0;bb<len;bb++){
		if(min>arr[bb]){
			min=arr[bb];
		}
	}
	return min;
}

double FindAverage(double * arr,int len){
	int cc;
	double average =0;
	for(cc=0;cc<len;cc++){
		average = average + arr[cc];
	}
	return (average)/len;
}

int main(){
    int maxval,minval,len,i,x,y;
    double tstart,tstop,tdiff;
    maxval=500;
    minval=-50;
	double * quickAverage =(double *) malloc(10*sizeof(double));
    double * quickMin =(double *) malloc(10*sizeof(double));
    double * quickMax =(double *) malloc(10*sizeof(double));
	double * arrayDiff =(double *) malloc(10*sizeof(double));
    int a;
    for(a=0,len=16000;a<10;len=len*1.2,a++){
    	   	
    	for(x=0;x<10;x++){
    		int *array = random_array(maxval,minval,len);
    		reverse_sort(array,len);
    		tstart = clockms();
    		quick_sort(array,len);
    		tstop = clockms();
    		tdiff = tstop-tstart;
    		arrayDiff[x] = tdiff;

    		free(array);
		
		if(x==9){
			//	quickMin[a] = FindMin(arrayDiff,10); //sorun find fonksiyonlarýnda
    		/*		if(a==0){
			
    			quickAverage[a] = FindAverage(arrayDiff,10);
    			}
    			else{
    				quickAverage[a] = ((FindAverage(arrayDiff,10))/(pow(1.2,a)*log(pow(1.2,a))));
				}*/
    		//	quickMax[a] = FindMax(arrayDiff,10);
    		//	printf("quickMin[%d] = %f\n",a,quickMin[a]);
    		quickAverage[a] = FindAverage(arrayDiff,10);
	  			printf("%f\n",quickAverage[a]);
	  		//	printf("quickMax[%d] = %f\n",a,quickMax[a]);	  			
			}	
		}
	
	
	}
}
