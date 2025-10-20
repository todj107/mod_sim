#include <sys/types.h>
#include <sys/wait.h>
#include <stdio.h>
#include <errno.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>

void run_programs(char *program, char *args[], pid_t *pids, char *vbl[], int n, int idx){
	// Create a child process for each command
	for (int i = 0; i < n; i++) {
		args[idx] = vbl[i];
		pid_t pid = fork();
		if (pid == -1) {
			perror("Fork");
			exit(EXIT_FAILURE);
		}
		else if (pid == 0) {
			if(execvp(program, args) < 0){
				perror("Execvp");
				exit(EXIT_FAILURE);
			}
		}
		else {
            // Parent
            pids[i] = pid;
        }
	}

	// Wait for children
	for (size_t i = 0; i < n; i++){
		int status;
		if (waitpid(pids[i], &status, 0) == -1) {
			perror("Waitpid");
		}
		else if (WIFEXITED(status)) {
			int exit_code = WEXITSTATUS(status);
			if (exit_code != 0) {
				fprintf(stderr, "Child %zu (%s): Returned failure!\n", i, vbl[i]);
			}
		}
	}
}


int main(void){
	char *program = "./sim";

	// TASK 1
	char *args_task1[] = {program, "N=64", "rho=0.5", "T=1.0", "alpha=0", NULL, "nblock=20", "read=0064_r0.500_T1.000_start", "run", NULL};
	char *dt[] = {"deltat=0.002", "deltat=0.004", "deltat=0.006", "deltat=0.008", "deltat=0.010", "deltat=0.012"};
	int n1 = 6; // Number of executions.
	int idx_1 = 5;

	// TASK 2
	char *args_task2[] = {program, "N=64", "rho=0.6", "T=1.0", NULL, "deltat=0.01", "nblock=10", "run", NULL};
	char *alpha[] = {"alpha=0.01", "alpha=0.1", "alpha=1"};
	int n2 = 3;
	int idx_2 = 4;


	int n = n2;
	pid_t *pids = malloc(n * sizeof(pid_t)); 
	if (pids == NULL){
		perror("Malloc");
		exit(EXIT_FAILURE);
	}

	// run_programs(program, args_task1, pids, dt, n1, idx_1);
	run_programs(program, args_task2, pids, alpha, n2, idx_2);
	

	free(pids);
	return 0;
}