# PMNS for efficiency and small memory cost

Contains a generator of low memory efficient PMNS and an efficient C codes generator.

The PMNS generator in the directory 'pmns_generator' is based on this work: PMNS for efficient arithmetic and small memory cost (to appear in Advances in IEEE Transactions on Emerging Topics in Computing, 2022).

A standard C codes generator is provied in the directory 'C_codes_generator'. Given a PMNS (obtained from the preceding generator), it generates C codes for all the main operations (including forward and backward conversion to PMNS, addition, multiplication).
