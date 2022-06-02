#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <time.h>
#include <gmp.h>
#include <openssl/bn.h>

//~ #include "intel_measurement_stuff.c"
#include "gmp_stuff.c"

#include "structs_data.h"
#include "add_mult_poly.c"
#include "useful_functs.c"


#define BILLION 1000000000L


//~ Compilation & execution command: gcc -Wall -O3 main.c -o main -lgmp -lcrypto && ./main

//~ Important : polynomials representations form is P(X) = a0 + ... + an.X^n = (a0, ..., an).



int main(void){
	
	srand(time(NULL));	
	
	unsigned long seed = time(NULL);
	gmp_randstate_t r;
	gmp_randinit_default(r);
	gmp_randseed_ui(r, seed);

	BN_CTX *ctx = BN_CTX_new();
	BN_CTX_start(ctx);

	BIGNUM *opA = BN_CTX_get(ctx);
	BIGNUM *opB = BN_CTX_get(ctx);
	BIGNUM *opP = BN_CTX_get(ctx);
	BIGNUM *opAA = BN_CTX_get(ctx);
	BIGNUM *opBB = BN_CTX_get(ctx);
	BN_MONT_CTX *mont_ctx = BN_MONT_CTX_new();

	struct timespec start0, end0;
	struct timespec start1, end1;
	struct timespec start2, end2;
	struct timespec start3, end3;
	struct timespec start4, end4;
	struct timespec start5, end5;
	struct timespec start6, end6;
	struct timespec start7, end7;
	struct timespec start8, end8;
	struct timespec start9, end9;
	struct timespec start10, end10;
	uint64_t diff0, diff1, diff2, diff3, diff4, diff5, diff6;
	uint64_t diff7, diff8, diff9, diff10;

	int i, nbiter, nb_limbs;
	mpz_t A, B, C, E, F, R2, G, H;
	mpz_inits (A, B, C, E, F, R2, G, H, NULL);
	
	mp_limb_t mip0;
	
	const mp_limb_t *p_limbs;
	mp_limb_t *r2_limbs, *scratch_limbs, *mip_limbs;
	mp_limb_t *a_limbs, *b_limbs, *am1_limbs, *am2_limbs, *bm1_limbs, *bm2_limbs;
	
	int64_t pa[NB_COEFF];
	int64_t pb[NB_COEFF];
	
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	init_data();
	
	mpz_urandomm(A, r, modul_p);
	mpz_urandomm(B, r, modul_p);
	mpz_set(E, A);
	
	nb_limbs = mpz_size (modul_p);
	
	p_limbs = mpz_limbs_read (modul_p);   
	
	binvert_limb (mip0, p_limbs[0]);
	mip0 = -mip0;
	
	mpz_setbit (R2, 2*nb_limbs*8*sizeof(mp_limb_t)); 
	mpz_mod(R2, R2, modul_p);
	r2_limbs = mpz_limbs_modify (R2, nb_limbs); 
	
	a_limbs = (mp_limb_t*) calloc (nb_limbs, sizeof(mp_limb_t));
	am1_limbs = (mp_limb_t*) calloc (nb_limbs, sizeof(mp_limb_t));
	am2_limbs = (mp_limb_t*) calloc (nb_limbs, sizeof(mp_limb_t));
	bm1_limbs = (mp_limb_t*) calloc (nb_limbs, sizeof(mp_limb_t));
	bm2_limbs = (mp_limb_t*) calloc (nb_limbs, sizeof(mp_limb_t));
	mip_limbs = (mp_limb_t*) calloc (nb_limbs, sizeof(mp_limb_t));
	scratch_limbs = (mp_limb_t*) calloc (2*nb_limbs, sizeof(mp_limb_t));

	b_limbs = mpz_limbs_modify (B, nb_limbs);
	copy_limbs(a_limbs, A, nb_limbs);
	copy_limbs(am1_limbs, A, nb_limbs);
	copy_limbs(am2_limbs, A, nb_limbs);
	copy_limbs(bm1_limbs, B, nb_limbs);
	copy_limbs(bm2_limbs, B, nb_limbs);
	
	mpn_binvert (mip_limbs, p_limbs, nb_limbs, scratch_limbs); //clean_limbs(scratch_limbs, 2*nb_limbs);
	
	BN_dec2bn(&opA, mpz_get_str (NULL, 10, A));
	BN_dec2bn(&opB, mpz_get_str (NULL, 10, B));
	BN_dec2bn(&opP, mpz_get_str (NULL, 10, modul_p));
	BN_copy(opAA, opA);
	BN_copy(opBB, opB);
	BN_MONT_CTX_set(mont_ctx, opP, ctx);
	
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID , &start3);
	from_int_to_amns(pa, A);
	from_int_to_amns(pb, B);
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID , &end3);
	diff3 = BILLION * (end3.tv_sec - start3.tv_sec) + (end3.tv_nsec - start3.tv_nsec);

	clock_gettime(CLOCK_PROCESS_CPUTIME_ID , &start5);
	BN_to_montgomery(opA, opA, mont_ctx, ctx);
	BN_to_montgomery(opB, opB, mont_ctx, ctx);
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID , &end5);
	diff5 = BILLION * (end5.tv_sec - start5.tv_sec) + (end5.tv_nsec - start5.tv_nsec);
	
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID , &start9);
	//~ conversion to Montgomery domain (mont par bloc)
	mpn_mont_mul_red_1(am1_limbs, am1_limbs, r2_limbs, p_limbs, mip0, nb_limbs);
	mpn_mont_mul_red_1(bm1_limbs, bm1_limbs, r2_limbs, p_limbs, mip0, nb_limbs);
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID , &end9);
	diff9 = BILLION * (end9.tv_sec - start9.tv_sec) + (end9.tv_nsec - start9.tv_nsec);
	
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID , &start10);
	//~ conversion to Montgomery domain (mont classique)
	mpn_mont_mul_red_n(am2_limbs, am2_limbs, r2_limbs, p_limbs, mip_limbs, nb_limbs);
	mpn_mont_mul_red_n(bm2_limbs, bm2_limbs, r2_limbs, p_limbs, mip_limbs, nb_limbs);
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID , &end10);
	diff10 = BILLION * (end10.tv_sec - start10.tv_sec) + (end10.tv_nsec - start10.tv_nsec);

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	nbiter = 1 << 2;
	
	
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID , &start1);
	for (i=0; i<nbiter; i++) {
		mpz_mul (A, A, B);
		mpz_mod (A, A, modul_p);
	}
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID , &end1);
	diff1 = BILLION * (end1.tv_sec - start1.tv_sec) + (end1.tv_nsec - start1.tv_nsec);
	
	
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID , &start0);
	for (i=0; i<nbiter; i++) {
		mpn_mod_mult(a_limbs, a_limbs, b_limbs, p_limbs, nb_limbs);
	}
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID , &end0);
	diff0 = BILLION * (end0.tv_sec - start0.tv_sec) + (end0.tv_nsec - start0.tv_nsec);

	
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID , &start7);
	for (i=0; i<nbiter; i++) {
		//~ Montgomery modular multiplication (mont par bloc)
		mpn_mont_mul_red_1(am1_limbs, am1_limbs, bm1_limbs, p_limbs, mip0, nb_limbs);
	}
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID , &end7);
	diff7 = BILLION * (end7.tv_sec - start7.tv_sec) + (end7.tv_nsec - start7.tv_nsec);
	
	
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID , &start8);
	for (i=0; i<nbiter; i++) {
		//~ Montgomery modular multiplication (mont classique)
		mpn_mont_mul_red_n(am2_limbs, am2_limbs, bm2_limbs, p_limbs, mip_limbs, nb_limbs);
	}
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID , &end8);
	diff8 = BILLION * (end8.tv_sec - start8.tv_sec) + (end8.tv_nsec - start8.tv_nsec);
	

	clock_gettime(CLOCK_PROCESS_CPUTIME_ID , &start2);
	for (i=0; i<nbiter; i++) {
		mult_mod_poly(pa, pa, pb);
	}
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID , &end2);
	diff2 = BILLION * (end2.tv_sec - start2.tv_sec) + (end2.tv_nsec - start2.tv_nsec);


	clock_gettime(CLOCK_PROCESS_CPUTIME_ID , &start4);
	for (i=0; i<nbiter; i++) {
		BN_mod_mul_montgomery(opA, opA, opB, mont_ctx, ctx);
	}
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID , &end4);
	diff4 = BILLION * (end4.tv_sec - start4.tv_sec) + (end4.tv_nsec - start4.tv_nsec);


	clock_gettime(CLOCK_PROCESS_CPUTIME_ID , &start6);
	for (i=0; i<nbiter; i++) {
		BN_mod_mul(opAA, opAA, opBB, opP, ctx);
	}
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID , &end6);
	diff6 = BILLION * (end6.tv_sec - start6.tv_sec) + (end6.tv_nsec - start6.tv_nsec);
	

	//~ should not modify the value which is represented
	exact_coeffs_reduction(pa, pa);

	clock_gettime(CLOCK_PROCESS_CPUTIME_ID , &start3);
	from_amns_to_int(C, pa);
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID , &end3);
	diff3 += BILLION * (end3.tv_sec - start3.tv_sec) + (end3.tv_nsec - start3.tv_nsec);

	clock_gettime(CLOCK_PROCESS_CPUTIME_ID , &start5);
	BN_from_montgomery(opA, opA, mont_ctx, ctx);
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID , &end5);
	diff5 += BILLION * (end5.tv_sec - start5.tv_sec) + (end5.tv_nsec - start5.tv_nsec);
	
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	from_limbs_to_mpz_t(F, a_limbs, nb_limbs);
	
	clean_limbs(scratch_limbs, 2*nb_limbs);
	for(i=0; i<nb_limbs; i++)
		scratch_limbs[i] = am1_limbs[i];
	
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID , &start9);
	mpn_redc_1 (am1_limbs, scratch_limbs, p_limbs, nb_limbs, mip0);
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID , &end9);
	diff9 += BILLION * (end9.tv_sec - start9.tv_sec) + (end9.tv_nsec - start9.tv_nsec);
	
	from_limbs_to_mpz_t(G, am1_limbs, nb_limbs);
	
	
	clean_limbs(scratch_limbs, 2*nb_limbs);
	for(i=0; i<nb_limbs; i++)
		scratch_limbs[i] = am2_limbs[i];
	
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID , &start10);
	mpn_redc_n (am2_limbs, scratch_limbs, p_limbs, nb_limbs, mip_limbs);
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID , &end10);
	diff10 += BILLION * (end10.tv_sec - start10.tv_sec) + (end10.tv_nsec - start10.tv_nsec);
	
	from_limbs_to_mpz_t(H, am2_limbs, nb_limbs);

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	printf("\nnbiter = %d\n\n", nbiter);
	gmp_printf("p       : %Zd\n\n", modul_p);
	gmp_printf("A       : %Zd\n", E);
	gmp_printf("B       : %Zd\n\n", B);
	gmp_printf("r_gmp   : %Zd\n", A);
	gmp_printf("r_lgmp  : %Zd\n", F);
	gmp_printf("r_mbgmp : %Zd\n", G);
	gmp_printf("r_mcgmp : %Zd\n", H);
	gmp_printf("r_amns  : %Zd\n", C);
	gmp_printf("r_ssld  : %s\n", BN_bn2dec(opAA));
	gmp_printf("r_sslm  : %s\n\n\n", BN_bn2dec(opA));
	
	//~ printf("Modular multiplication mean timings (in nanoseconds)\n\n");
	printf("      ********************************************************\n");
	printf("      * Modular multiplication mean timings (in nanoseconds) *\n");
	printf("      ********************************************************\n\n");
	printf("gnu mp		= %lu\n", (diff1/nbiter));
	printf("low gnu mp	= %lu\n", (diff0/nbiter));
	printf("bmont gnu mp	= %lu\n", (diff7/nbiter));
	printf("cmont gnu mp	= %lu\n", (diff8/nbiter));
	printf("default openssl	= %lu\n", (diff6/nbiter));
	printf("mont openssl	= %lu\n", (diff4/nbiter));
	printf("amns		= %lu\n\n", (diff2/nbiter));
	
	printf("conv time for amns		= %lu\n", diff3);
	printf("conv time for openssl mont	= %lu\n", diff5);
	printf("conv time for gmp bloc mont	= %lu\n", diff9);
	printf("conv time for gmp clas mont	= %lu\n\n", diff10);


	mpz_clears (A, B, C, E, F, R2, G, H, NULL);
	gmp_randclear(r);

	BN_MONT_CTX_free(mont_ctx);
	BN_CTX_end(ctx);
	BN_CTX_free(ctx);
	
	free(a_limbs);
	free(am1_limbs);
	free(am2_limbs);
	free(bm1_limbs);
	free(bm2_limbs);
	free(mip_limbs);
	free(scratch_limbs);
	
	free_data();
	return 0;
}

















