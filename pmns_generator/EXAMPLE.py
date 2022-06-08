load("pmns_generator.sage")

roots_max_duration_checks = 60  # seconds

bSize = 256
p = random_prime(2**bSize, lbound=2**(bSize-1))
# or p = Integer(...) to give a specific value for p

nb_free_add = 1
phi_log2 = 64 
for_avx = False

lambda_max = 2

n = floor(p.nbits()/phi_log2) + 1


pmns_list = build_pmns_candidates_for_n(p, n, phi_log2, nb_free_add, lambda_max, roots_max_duration_checks, for_avx)


if len(pmns_list) != 0 : print(pmns_list[0])


for pmns in pmns_list:
	print(pmns)
	print()



#~ NOTE: data struct: [nb_max_add, n, E, rhoUp_log2, phi_log2, gmm, M]
