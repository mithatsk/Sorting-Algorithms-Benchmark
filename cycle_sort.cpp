#include<stdio.h>
#include <time.h>
#include <stdlib.h> 
#include <sys/time.h>
#include<math.h>
void swap(int *,int *);
void gnome_sort(int *,int);
int * random_array(int,int,int);
static double clockms (void);
double FindMax(double *,int);
double FindMin(double *,int);
double FindAverage(double *,int);

void cycle_sort(int * list, size_t l_len);
void show_array(int * array, size_t a_len);
 

void cycle_sort(int * list, size_t l_len)
{
  
 
  /* Loop through the array to find cycles to rotate. */
  for (int cycleStart = 0; cycleStart < l_len - 1; ++cycleStart)
  {
    int item = list[cycleStart];
    int swap_tmp;
 
    /* Find where to put the item. */
    int pos = cycleStart;
    for (int i = cycleStart + 1; i < l_len; ++i)
    {
      if (list[i] < item)
      {
        ++pos;
      }
    }
 
    /* If the item is already there, this is not a cycle. */
    if (pos == cycleStart)
    {
      continue;
    }
 
    /* Otherwise, put the item there or right after any duplicates. */
    while (item == list[pos])
    {
      ++pos;
    }
    swap_tmp = list[pos];
    list[pos] = item;
    item = swap_tmp;
  
 
    /* Rotate the rest of the cycle. */
    while (pos != cycleStart)
    {
      /* Find where to put the item. */
      pos = cycleStart;
      for (int i = cycleStart + 1; i < l_len; ++i)
      {
        if (list[i] < item)
        {
          ++pos;
        }
      }
 
      /* Put the item there or right after any duplicates. */
      while (item == list[pos])
      {
        ++pos;
      }
      swap_tmp = list[pos];
      list[pos] = item;
      item = swap_tmp;
  
    }
  }

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



void show_array(int * array, size_t a_len)
{
  for (int ix = 0; ix < a_len; ++ix)
  {
    printf("%d ", array[ix]);
  }
  putchar('\n');
 
  return;
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
	double * cycleAverage =(double *) malloc(10*sizeof(double));
    double * cycleMin =(double *) malloc(10*sizeof(double));
    double * cycleMax =(double *) malloc(10*sizeof(double));
	double * arrayDiff =(double *) malloc(10*sizeof(double));
    int a;
    for(a=0,len=16000;a<10;len=len*1.2,a++){
    	   	
    	for(x=0;x<10;x++){
    		int *array = random_array(maxval,minval,len);
    		reverse_sort(array,len);
    		tstart = clockms();
    		cycle_sort(array,len);
    		tstop = clockms();
    		tdiff = tstop-tstart;
    		arrayDiff[x] = tdiff;

    		free(array);
		
		if(x==9){
			//	cycleMin[a] = FindMin(arrayDiff,10); //sorun find fonksiyonlarýnda
			cycleAverage[a] = FindAverage(arrayDiff,10);
    		//	cycleAverage[a] = (FindAverage(arrayDiff,10))/pow(1.2,2*a);
    		//	cycleMax[a] = FindMax(arrayDiff,10);
    		//	printf("cycleMin[%d] = %f\n",a,cycleMin[a]);
	  			printf("%f\n",cycleAverage[a]);
	  		//	printf("cycleMax[%d] = %f\n",a,cycleMax[a]);	  			
			}	
		}
	
	
	}
}
