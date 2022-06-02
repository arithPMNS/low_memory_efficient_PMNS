#include <unistd.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <time.h>
#include <gmp.h>
#include <openssl/bn.h>

#include "gmp_stuff.c"
#include "intel_measurement_stuff.c"

#include "structs_data.h"
#include "add_mult_poly.c"
#include "useful_functs.c"


//~ Compilation & execution command: gcc -Wall -O3 main.c -o main -lgmp -lcrypto && ./main

//~ Important: polynomials representations form is P(X) = a0 + ... + an.X^n = (a0, ..., an).


int main(void){
	
	srand(time(NULL));
	
	unsigned long seed = time(NULL);
	gmp_randstate_t r;
	gmp_randinit_default(r);
	gmp_randseed_ui(r, seed);
	
	unsigned long long timermin, timermax, t1, t2, diff_t;
	unsigned long long meanTimer1_min=0, meanTimer2_min=0, meanTimer3_min=0, meanTimer4_min=0, meanTimer5_min=0, meanTimer6_min=0;
	unsigned long long meanTimer1_max=0, meanTimer2_max=0, meanTimer3_max=0, meanTimer4_max=0, meanTimer5_max=0, meanTimer6_max=0;
	unsigned long long *statTimer1, *statTimer2, *statTimer3, *statTimer4, *statTimer5, *statTimer6;
	uint64_t cycles1[NTEST*NSAMPLES]={0}, cycles2[NTEST*NSAMPLES]={0}, cycles3[NTEST*NSAMPLES]={0}, cycles4[NTEST*NSAMPLES]={0}, cycles5[NTEST*NSAMPLES]={0}, cycles6[NTEST*NSAMPLES]={0};
	
	//~ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	int nb_limbs;
	mp_limb_t mip0;
	
	mpz_t R2, A[NTEST], B[NTEST], C[NTEST];
	const mp_limb_t *p_limbs, *a_limbs[NTEST], *b_limbs[NTEST];
	mp_limb_t *am_limbs[NTEST], *bm_limbs[NTEST], *rm_limbs[NTEST], *r_limbs[NTEST], *r2_limbs;
	
	BN_CTX *ctx = BN_CTX_new();
	BN_MONT_CTX *mont_ctx = BN_MONT_CTX_new();
	BIGNUM *opP, *opA[NTEST], *opB[NTEST], *opC[NTEST], *opAA[NTEST], *opBB[NTEST], *opCC[NTEST];
	
	int64_t pa[NTEST][NB_COEFF];
	int64_t pb[NTEST][NB_COEFF];
	int64_t pc[NTEST][NB_COEFF];
	
	//~ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	mpz_init (R2);
	init_data(); // <---- for AMNS and the module (modul_p)
	
	nb_limbs = mpz_size (modul_p);
	
	p_limbs = mpz_limbs_read (modul_p);
	
	binvert_limb (mip0, p_limbs[0]);
	mip0 = -mip0;
	
	mpz_setbit (R2, 2*nb_limbs*8*sizeof(mp_limb_t)); 
	mpz_mod(R2, R2, modul_p);
	r2_limbs = mpz_limbs_modify (R2, nb_limbs); 
	
	BN_CTX_start(ctx);
	opP = BN_CTX_get(ctx);
	BN_dec2bn(&opP, mpz_get_str (NULL, 10, modul_p));
	BN_MONT_CTX_set(mont_ctx, opP, ctx);
	
	for (int i=0; i<NTEST; i++){
		mpz_init (A[i]);
		mpz_init (B[i]);
		mpz_init (C[i]);
		
		opA[i] = BN_CTX_get(ctx);
		opB[i] = BN_CTX_get(ctx);
		opC[i] = BN_CTX_get(ctx);
		
		opAA[i] = BN_CTX_get(ctx);
		opBB[i] = BN_CTX_get(ctx);
		opCC[i] = BN_CTX_get(ctx);
		
		r_limbs[i] = (mp_limb_t*) calloc (nb_limbs, sizeof(mp_limb_t));
		rm_limbs[i] = (mp_limb_t*) calloc (nb_limbs, sizeof(mp_limb_t));
		am_limbs[i] = (mp_limb_t*) calloc (nb_limbs, sizeof(mp_limb_t));
		bm_limbs[i] = (mp_limb_t*) calloc (nb_limbs, sizeof(mp_limb_t));
	}
	
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	////////////////////////// cache memory heating ////////////////////
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	for(int i=0;i<NTEST;i++){
		mpz_urandomm(A[i], r, modul_p);
		mpz_urandomm(B[i], r, modul_p);
		
		a_limbs[i] = mpz_limbs_read (A[i]);
		b_limbs[i] = mpz_limbs_read (B[i]);
		
		copy_limbs(am_limbs[i], A[i], nb_limbs);
		copy_limbs(bm_limbs[i], B[i], nb_limbs);
		mpn_mont_mul_red_1(am_limbs[i], am_limbs[i], r2_limbs, p_limbs, mip0, nb_limbs);
		mpn_mont_mul_red_1(bm_limbs[i], bm_limbs[i], r2_limbs, p_limbs, mip0, nb_limbs);
		
		from_int_to_amns(pa[i], A[i]);
		from_int_to_amns(pb[i], B[i]);
		
		BN_dec2bn(&(opA[i]), mpz_get_str (NULL, 10, A[i]));
		BN_dec2bn(&(opB[i]), mpz_get_str (NULL, 10, B[i]));
		
		BN_copy(opAA[i], opA[i]);
		BN_copy(opBB[i], opB[i]);
		BN_to_montgomery(opAA[i], opAA[i], mont_ctx, ctx);
		BN_to_montgomery(opBB[i], opBB[i], mont_ctx, ctx);
	}
	
	for(int i=0;i<NTEST;i++){
		BN_mod_mul(opC[i], opA[i], opB[i], opP, ctx);
		
		BN_mod_mul_montgomery(opCC[i], opAA[i], opBB[i], mont_ctx, ctx);
		
		mult_mod_poly(pc[i], pa[i], pb[i]);
		
		mpz_mul (C[i], A[i], B[i]);
		mpz_mod (C[i], C[i], modul_p);
		
		mpn_mod_mult(r_limbs[i], a_limbs[i], b_limbs[i], p_limbs, nb_limbs);
		
		mpn_mont_mul_red_1(rm_limbs[i], am_limbs[i], bm_limbs[i], p_limbs, mip0, nb_limbs);
	}
	
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	////////////////////////// timing ////////////////////////////////
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	for(int i=0;i<NSAMPLES;i++){
		
		// Génération d'un jeu de paramètres aléatoires
		for(int j=0;j<NTEST;j++){
			mpz_urandomm(A[j], r, modul_p);
			mpz_urandomm(B[j], r, modul_p);
			
			a_limbs[j] = mpz_limbs_read (A[j]);
			b_limbs[j] = mpz_limbs_read (B[j]);
			
			copy_limbs(am_limbs[j], A[j], nb_limbs);
			copy_limbs(bm_limbs[j], B[j], nb_limbs);
			mpn_mont_mul_red_1(am_limbs[j], am_limbs[j], r2_limbs, p_limbs, mip0, nb_limbs);
			mpn_mont_mul_red_1(bm_limbs[j], bm_limbs[j], r2_limbs, p_limbs, mip0, nb_limbs);
			
			from_int_to_amns(pa[j], A[j]);
			from_int_to_amns(pb[j], B[j]);
			
			BN_dec2bn(&(opA[j]), mpz_get_str (NULL, 10, A[j]));
			BN_dec2bn(&(opB[j]), mpz_get_str (NULL, 10, B[j]));
			
			BN_copy(opAA[j], opA[j]);
			BN_copy(opBB[j], opB[j]);
			BN_to_montgomery(opAA[j], opAA[j], mont_ctx, ctx);
			BN_to_montgomery(opBB[j], opBB[j], mont_ctx, ctx);
		}
		
		//~~~~~~~~~~~~~~~~~~~~~~~~~~ for OpenSSL ~~~~~~~~~~~~~~~~~~~~~~~
		
		timermin = (unsigned long long int)0x1<<63;
		timermax = 0;
		for(int j=0;j<NTEST;j++){
			t1 = cpucyclesStart();
			
			BN_mod_mul(opC[j], opA[j], opB[j], opP, ctx); // Appel de la fonction à mesurer, avec le jeu de paramètres précédant
			
			t2 = cpucyclesStop();
			if (t2 < t1){
				diff_t = 18446744073709551615ULL-t1;
				diff_t = t2+diff_t+1;
			}
			else
				diff_t = t2-t1;
			
			if(timermin > diff_t) 
				timermin = diff_t;
			else if(timermax < diff_t) 
				timermax = diff_t;
			
			cycles1[i*NTEST+j]=diff_t;
		}

		meanTimer1_min += timermin;
		meanTimer1_max += timermax;
		
		
		//~~~~~~~~~~~~~~~~~~~~~~~~~~ for OpenSSL (MONTGOMERY) ~~~~~~~~~~
		
		timermin = (unsigned long long int)0x1<<63;
		timermax = 0;
		for(int j=0;j<NTEST;j++){
			t1 = cpucyclesStart();
			
			BN_mod_mul_montgomery(opCC[j], opAA[j], opBB[j], mont_ctx, ctx); // Appel de la fonction à mesurer, avec le jeu de paramètres précédant
			
			t2 = cpucyclesStop();
			if (t2 < t1){
				diff_t = 18446744073709551615ULL-t1;
				diff_t = t2+diff_t+1;
			}
			else
				diff_t = t2-t1;
			
			if(timermin > diff_t) 
				timermin = diff_t;
			else if(timermax < diff_t) 
				timermax = diff_t;
			
			cycles2[i*NTEST+j]=diff_t;
		}

		meanTimer2_min += timermin;
		meanTimer2_max += timermax;
		
		
		//~~~~~~~~~~~~~~~~~~~~~~~~~~ for AMNS ~~~~~~~~~~~~~~~~~~~~~~~~~~
		
		timermin = (unsigned long long int)0x1<<63;
		timermax = 0;
		for(int j=0;j<NTEST;j++){
			t1 = cpucyclesStart();
			
			mult_mod_poly(pc[j], pa[j], pb[j]); // Appel de la fonction à mesurer, avec le jeu de paramètres précédant
			
			t2 = cpucyclesStop();
			if (t2 < t1){
				diff_t = 18446744073709551615ULL-t1;
				diff_t = t2+diff_t+1;
			}
			else
				diff_t = t2-t1;
			
			if(timermin > diff_t) 
				timermin = diff_t;
			else if(timermax < diff_t) 
				timermax = diff_t;
			
			cycles3[i*NTEST+j]=diff_t;
		}

		meanTimer3_min += timermin;
		meanTimer3_max += timermax;
		
		//~~~~~~~~~~~~~~~~~~~~~~~~~~ for GNU MP ~~~~~~~~~~~~~~~~~~~~~~~~
		
		timermin = (unsigned long long int)0x1<<63;
		timermax = 0;
		for(int j=0;j<NTEST;j++){
			t1 = cpucyclesStart();
			
			// Appel de la fonction à mesurer, avec le jeu de paramètres précédant
			mpz_mul (C[j], A[j], B[j]);
			mpz_mod (C[j], C[j], modul_p);
			
			t2 = cpucyclesStop();
			if (t2 < t1){
				diff_t = 18446744073709551615ULL-t1;
				diff_t = t2+diff_t+1;
			}
			else
				diff_t = t2-t1;
			
			if(timermin > diff_t) 
				timermin = diff_t;
			else if(timermax < diff_t) 
				timermax = diff_t;
			
			cycles4[i*NTEST+j]=diff_t;
		}

		meanTimer4_min += timermin;
		meanTimer4_max += timermax;
		
		//~~~~~~~~~~~~~~~~~~~~~~~~~~ for GNU MP (low level) ~~~~~~~~~~~~
		
		timermin = (unsigned long long int)0x1<<63;
		timermax = 0;
		for(int j=0;j<NTEST;j++){
			t1 = cpucyclesStart();
			
			// Appel de la fonction à mesurer, avec le jeu de paramètres précédant
			mpn_mod_mult(r_limbs[j], a_limbs[j], b_limbs[j], p_limbs, nb_limbs);
			
			t2 = cpucyclesStop();
			if (t2 < t1){
				diff_t = 18446744073709551615ULL-t1;
				diff_t = t2+diff_t+1;
			}
			else
				diff_t = t2-t1;
			
			if(timermin > diff_t) 
				timermin = diff_t;
			else if(timermax < diff_t) 
				timermax = diff_t;
			
			cycles5[i*NTEST+j]=diff_t;
		}

		meanTimer5_min += timermin;
		meanTimer5_max += timermax;
		
		//~~~~~~~~~~~~~~~~~~~~~~~~~~ for GNU MP (MONTGOMERY) ~~~~~~~~~~~
		
		timermin = (unsigned long long int)0x1<<63;
		timermax = 0;
		for(int j=0;j<NTEST;j++){
			t1 = cpucyclesStart();
			
			// Appel de la fonction à mesurer, avec le jeu de paramètres précédant
			mpn_mont_mul_red_1(rm_limbs[j], am_limbs[j], bm_limbs[j], p_limbs, mip0, nb_limbs);
			
			t2 = cpucyclesStop();
			if (t2 < t1){
				diff_t = 18446744073709551615ULL-t1;
				diff_t = t2+diff_t+1;
			}
			else
				diff_t = t2-t1;
			
			if(timermin > diff_t) 
				timermin = diff_t;
			else if(timermax < diff_t) 
				timermax = diff_t;
			
			cycles6[i*NTEST+j]=diff_t;
		}

		meanTimer6_min += timermin;
		meanTimer6_max += timermax;
		
		
	}
	statTimer1 = quartiles(cycles1, NTEST*NSAMPLES);
	statTimer2 = quartiles(cycles2, NTEST*NSAMPLES);
	statTimer3 = quartiles(cycles3, NTEST*NSAMPLES);
	statTimer4 = quartiles(cycles4, NTEST*NSAMPLES);
	statTimer5 = quartiles(cycles5, NTEST*NSAMPLES);
	statTimer6 = quartiles(cycles6, NTEST*NSAMPLES);

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	////////////////////////// results ////////////////////////////////
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	printf("\n      ***********************************************\n");
	printf("      * Modular multiplication results (CPU cycles) *\n");
	printf("      ***********************************************\n\n");
	printf("                ---------------------------------------\n");
	printf("               |   min |   max |    Q1 |    Q2 |    Q3 |\n");
	printf(" ------------------------------------------------------\n");
	printf("| amns         | %5lld | %5lld | %5lld | %5lld | %5lld | \n", meanTimer3_min/NSAMPLES, meanTimer3_max/NSAMPLES, statTimer3[0], statTimer3[1], statTimer3[2]);
	printf(" ------------------------------------------------------\n");
	printf("| gnu_mp_mont  | %5lld | %5lld | %5lld | %5lld | %5lld |\n", meanTimer6_min/NSAMPLES, meanTimer6_max/NSAMPLES, statTimer6[0], statTimer6[1], statTimer6[2]);
	printf(" ------------------------------------------------------\n");
	printf("| openssl_mont | %5lld | %5lld | %5lld | %5lld | %5lld | \n", meanTimer2_min/NSAMPLES, meanTimer2_max/NSAMPLES, statTimer2[0], statTimer2[1], statTimer2[2]);
	printf(" ------------------------------------------------------\n");
	printf("| gnu_mp_low   | %5lld | %5lld | %5lld | %5lld | %5lld | \n", meanTimer5_min/NSAMPLES, meanTimer5_max/NSAMPLES, statTimer5[0], statTimer5[1], statTimer5[2]);
	printf(" ------------------------------------------------------\n");
	printf("| gnu_mp       | %5lld | %5lld | %5lld | %5lld | %5lld | \n", meanTimer4_min/NSAMPLES, meanTimer4_max/NSAMPLES, statTimer4[0], statTimer4[1], statTimer4[2]);
	printf(" ------------------------------------------------------\n");
	printf("| openssl      | %5lld | %5lld | %5lld | %5lld | %5lld | \n", meanTimer1_min/NSAMPLES, meanTimer1_max/NSAMPLES, statTimer1[0], statTimer1[1], statTimer1[2]);
	printf(" ------------------------------------------------------\n\n");
	
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	////////////////////////// cleaning ////////////////////////////////
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	mpz_clear (R2);
	for (int i=0; i<NTEST; i++){
		mpz_clear (A[i]);
		mpz_clear (B[i]);
		mpz_clear (C[i]);
		
		free(r_limbs[i]);
		free(rm_limbs[i]);
		free(am_limbs[i]);
		free(bm_limbs[i]);
	}
	gmp_randclear(r);
	
	BN_MONT_CTX_free(mont_ctx);
	BN_CTX_end(ctx);
	BN_CTX_free(ctx);
	
	free_data();
	
	return 0;
}









