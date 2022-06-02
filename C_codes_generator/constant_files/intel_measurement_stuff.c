#define NTEST 1001
#define NSAMPLES 51
// NTEST*NSAMPLES must be odd
// it's easier to compute median value


/**** Measurements procedures according to INTEL white paper

  "How to benchmark code execution times on INTEL IA-32 and IA-64"

 *****/

void quicksort(uint64_t *tab,int first,int last);
unsigned long long *quartiles(uint64_t *tab, uint64_t size);

void quicksort(uint64_t *tab,int first,int last){

	int i, j, pivot;
	uint64_t temp;

	if(first<last){
		pivot=first;
		i=first;
		j=last;

		while(i<j){
			while(tab[i]<=tab[pivot]&&i<last)
				i++;
			while(tab[j]>tab[pivot])
				j--;
			if(i<j){
				temp=tab[i];
				tab[i]=tab[j];
				tab[j]=temp;
			}
		}
		temp=tab[pivot];
		tab[pivot]=tab[j];
		tab[j]=temp;
		quicksort(tab,first,j-1);
		quicksort(tab,j+1,last);
	}
}

unsigned long long *quartiles(uint64_t *tab, uint64_t size)
{
	unsigned long long *result ;
	uint64_t aux ;

	result = malloc(3*sizeof(unsigned long long));
	quicksort(tab,0,size-1);
	aux = size >> 2;
	if (size % 4) aux++;
	// Q1
	result[0] = tab[aux-1];
	// Mediane
	// size is odd hence it's easy
	result[1]  = tab[(size+1)/2 - 1];
	// Q3
	aux = (3*size) >> 2;
	if ((3*size) % 4) aux++;
	result[2]  = tab[aux - 1];

	return result;
}

inline static uint64_t cpucyclesStart (void) {

	unsigned hi, lo;
	__asm__ __volatile__ (	"CPUID\n\t"
			"RDTSC\n\t"
			"mov %%edx, %0\n\t"
			"mov %%eax, %1\n\t"
			: "=r" (hi), "=r" (lo)
			:
			: "%rax", "%rbx", "%rcx", "%rdx");

	return ((uint64_t)lo)^(((uint64_t)hi)<<32);


}

inline static uint64_t cpucyclesStop (void) {

	unsigned hi, lo;
	__asm__ __volatile__(	"RDTSCP\n\t"
			"mov %%edx, %0\n\t"
			"mov %%eax, %1\n\t"
			"CPUID\n\t"
			: "=r" (hi), "=r" (lo)
			:
			: "%rax", "%rbx", "%rcx", "%rdx");

	return ((uint64_t)lo)^(((uint64_t)hi)<<32);


}

