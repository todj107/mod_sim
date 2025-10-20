#include <sys/types.h>
#include <sys/wait.h>
#include <stdio.h>
#include <errno.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>

void run_programs(char *program, char *args[], pid_t *pids, char *vbl[], int n){
	// Create a child process for each command
	for (int i = 0; i < n; i++) {
		args[5] = vbl[i];
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

	// TASK 5
	char *args_brown[] = {"./sim", "N=64", "rho=0.3", "T=1.0", "alpha=10", NULL, "nblock=100", "run", NULL};
	char *dt[] = {"deltat=0.001", "deltat=0.002", "deltat=0.003", "deltat=0.004"};
	int n = 4; // Number of executions.

	pid_t *pids = malloc((n) * sizeof(pid_t)); 
	if (pids == NULL){
		perror("Malloc");
		exit(EXIT_FAILURE);
	}

	run_programs(program, args_brown, pids, dt, n);
	

	free(pids);
	return 0;
}