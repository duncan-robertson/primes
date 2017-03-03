/* factorc
 * Prime Factoring
 * Duncan Robertson
 */

#include<stdio.h>
#include<stdlib.h>
#include<sys/resource.h>
#include<gmp.h>
#include<mpi.h>

MPI_Status status;
MPI_Request *reqs;
mpz_t prime1, prime2, n, tmp1, tmp2, quotient, rem, sz;
mpz_t *section;
unsigned int found=0, i, fail=0, bsize=0, startup=1;
int cmp1 = 0, my_rank, p, success=0;
char *buffer;

double getTime() {
    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);
    struct timeval time;
    time = usage.ru_utime;
    return time.tv_sec+time.tv_usec/1e6;
}

void factorCheck() {
    for(i=0; i<p; i++) {
        //Start new section
        mpz_set(prime1, section[i]);
        mpz_sub_ui(prime1, prime1, 1);

        mpz_nextprime(prime1, prime1);
        if(found)
            return;

        //Only select primes from your section, unless you are working on the final seciton
        //If working on final section ensure prime1 squared is less than n
        if(my_rank == p-1 && i == p-1) {
            mpz_mul(tmp1, prime1, prime1);
            cmp1 = mpz_cmp(tmp1, n);
        }
        else {
            mpz_sub(tmp1, prime1, section[i]);
            cmp1 = mpz_cmp(tmp1, sz);
        }

        while(cmp1 < 0) {
            //Divide n by prime1
            mpz_tdiv_qr(quotient, rem, n, prime1);

            //If n is cleanly divisible check if the quotient is prime
            if(mpz_cmp_ui(rem, 0) == 0) {
                switch(mpz_probab_prime_p(quotient, 15)) {
                    case(2):
                        mpz_set(prime2, quotient);
                        success = 1;
                        break;
                    case(1):
                        mpz_set(prime2, quotient);
                        mpz_sub_ui(prime2, prime2, 1);
                        mpz_nextprime(prime2, prime2);
                        if(mpz_cmp(prime2, quotient) == 0)
                            success = 1;
                        break;
                }

                //If it is prime, signal other processes
                if(success) {
                    if(my_rank < (p-1) && startup) {
                        bsize = 0;

                        for(i=1; i<=(p-1)-my_rank; i++) {
                            MPI_Send(&bsize, sizeof(unsigned int), MPI_CHAR, (my_rank+i), 1, MPI_COMM_WORLD);
                        }
                    }

                    MPI_Request r[p-1];

                    for(i=0; i < (p-1); i++) {
                        if(i < my_rank)
                            MPI_Isend(&success, 1, MPI_INT, i, 3, MPI_COMM_WORLD, &r[i]);
                        else 
                            MPI_Isend(&success, 1, MPI_INT, (i+1), 3, MPI_COMM_WORLD, &r[i]);
                    }

                    //Print primes, then leave function
                    gmp_printf("Prime factors found: %Zd %Zd\n", prime1, prime2);
                    return;
                }

            }

            mpz_nextprime(prime1, prime1);

            //If another process found the prime factors leave function
            if(found)
                return;

            if(my_rank == p-1 && i == p-1) {
                mpz_mul(tmp1, prime1, prime1);
                cmp1 = mpz_cmp(tmp1, n);
            }
            else {
                mpz_sub(tmp1, prime1, section[i]);
                cmp1 = mpz_cmp(tmp1, sz);
            }
        }
        //Section complete move to next section
    }
    //All sections complete, exit function
}

int main(int argc, char **argv) {
    //Initialize Variables
    char out[20];
    char filename[100];
    double time;
    FILE *f;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    reqs = malloc(sizeof(MPI_Request)*(p-1));
    section = malloc(sizeof(mpz_t)*p);

    mpz_init(prime1), mpz_init(prime2), mpz_init(n), mpz_init(tmp1), mpz_init(tmp2), mpz_init(quotient), mpz_init(rem), mpz_init(sz);

    for(i=0; i<p; i++) {
        mpz_init(section[i]);
    }

    sprintf(filename, "time_%s", argv[1]);

    //Validate argument
    if(my_rank == 0) {
        if(argc<2 || argc>3) {
            printf("Please supply one number to search for prime factors\n");
            fail = 1;
        }

        mpz_set_str(n, argv[1], 10);

        if(mpz_cmp_ui(n, 0) <= 0) {
            printf("Invalid number supplied.\nNon-numerical characters are illegal\n0 or negative numbers are invalid numbers\n");
            fail = 1;
        }
        f = fopen(filename, "w");
        if(f == NULL) {
            printf("There was a problem opening the file that would be used for output\n");
            fail = 1;
        }
        fclose(f);
    }

    //Broadcast failure status
    MPI_Bcast(&fail, 1, MPI_INT, 0, MPI_COMM_WORLD);

    //If something failed quit early
    if(fail){
        MPI_Finalize();
        exit(0);
    }

    if(my_rank != 0)
        mpz_set_str(n, argv[1], 10);

    //Computation started set starting time
    time = getTime();

    if(p == 1) {
        mpz_set_ui(prime1, 2);
        mpz_mul(tmp1, prime1, prime1);
        cmp1 = mpz_cmp(tmp1, n);

        //If prime1 squared is greater than n we checked all possible prime factors
        while(cmp1 < 0) {
            //Divide n by current prime
            mpz_tdiv_qr(quotient, rem, n, prime1);

            //If it divides evenly check if the quotient is prime
            if(mpz_cmp_ui(rem, 0) == 0) {
                switch(mpz_probab_prime_p(quotient, 15)){
                    case(2):
                        mpz_set(prime2, quotient);
                        found = 1;
                        break;
                    case(1):
                        mpz_set(prime2, quotient);
                        mpz_sub_ui(prime2, prime2, 1);
                        mpz_nextprime(prime2, prime2);
                        if(mpz_cmp(prime2, quotient) == 0)
                            found = 1;
                        break;
                }
            }

            //If the quotient was prime prime factors were found, halt computation
            if(found)
                break;

            //If the quotient wasn't prime check next possible prime
            mpz_nextprime(prime1, prime1);
            mpz_mul(tmp1, prime1, prime1);
            cmp1 = mpz_cmp(tmp1, n);
        }

        time = getTime()-time;

        sprintf(out, "%d\t% .2e\n", my_rank, time);

        do {
            f = fopen(filename, "a");
        } while (f == NULL);

        fputs(out, f);
        fclose(f);

        if(found)
            gmp_printf("Prime factors found %Zd %Zd\n", prime1, prime2);
        else
            printf("Prime factors could not be found\n");
    }
    else {
        if(my_rank == 0) {
            //Find the approximate square root of n
            mpz_sqrt(sz, n);

            //Divide that square root by number of processors squared
            mpz_cdiv_q_ui(sz, sz, p*p);

            //Send that size to other processes
            bsize = mpz_sizeinbase(sz, 2) + 2;
            buffer = malloc(sizeof(char)*bsize);
            mpz_get_str(buffer, 2, sz);
            for(i=1; i<p; i++) {
                MPI_Send(&bsize, sizeof(unsigned int), MPI_CHAR, i, 1, MPI_COMM_WORLD);
                MPI_Send(buffer, (sizeof(char)*bsize), MPI_CHAR, i, 2, MPI_COMM_WORLD);
            }
            free(buffer);
        }
        else {
            //Receive the size of the segment
            MPI_Recv(&bsize, sizeof(unsigned int), MPI_CHAR, 0, 1, MPI_COMM_WORLD, &status);

            buffer = malloc(sizeof(char)*bsize);
            MPI_Recv(buffer, (sizeof(char)*bsize), MPI_CHAR, 0, 2, MPI_COMM_WORLD, &status);

            mpz_set_str(sz, buffer, 2);
            free(buffer);
        }

        //Section i is equal to the start of a computation section
        for(i=0; i<p; i++) {
            mpz_mul_ui(tmp1, sz, i*p);
            mpz_mul_ui(tmp2, sz, my_rank);
            mpz_add(section[i], tmp1, tmp2);
        }

        //Set up nonblocking receive for success signal
        for(i=0; i < (p-1); i++) {
            if(i < my_rank)
                MPI_Irecv(&found, 1, MPI_INT, i, 3, MPI_COMM_WORLD, &reqs[i]);
            else 
                MPI_Irecv(&found, 1, MPI_INT, (i+1), 3, MPI_COMM_WORLD, &reqs[i]);
        }

        //Run primary computation
        factorCheck();

        //Computation done, write execution time to file
        time = getTime()-time;

        sprintf(out, "%d\t% .2e\n", my_rank, time);

        do {
            f = fopen(filename, "a");
        } while (f == NULL);

        fputs(out, f);
        fclose(f);
    }

    //Cleanup
    MPI_Finalize();
    
    mpz_clear(prime1), mpz_clear(prime2), mpz_clear(n), mpz_clear(tmp1), mpz_clear(tmp2), mpz_clear(quotient), mpz_clear(rem), mpz_clear(sz);
    
    for(i=0; i<p; i++) {
        mpz_clear(section[i]);
    }

    free(section);
    free(reqs);

    return 0;
}
