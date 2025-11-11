/* bank_simulator.c
   Clash of Coders - Bank Queue (Poisson) Simulator
   - Simulates customer arrivals (Poisson) for an 8-hour day (480 minutes)
   - Uses a linked-list queue implemented with structs, malloc/free
   - Simulates one or more tellers (service time 2-3 minutes)
   - Collects wait times dynamically and computes mean, median, mode, std-dev, max
   Compile: gcc bank_simulator.c -o bank_simulator -lm
*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>

/* -------------------------
   Structures and queue ops
   -------------------------*/

typedef struct Customer {
    int arrival_time;       // minute customer arrived (0..479)
    int service_start;      // minute service started
    int service_duration;   // how many minutes teller will take
    struct Customer *next;
} Customer;

typedef struct Queue {
    Customer *head;
    Customer *tail;
    int size;
} Queue;

void init_queue(Queue *q) {
    q->head = q->tail = NULL;
    q->size = 0;
}

void enqueue(Queue *q, Customer *c) {
    c->next = NULL;
    if (q->tail == NULL) {
        q->head = q->tail = c;
    } else {
        q->tail->next = c;
        q->tail = c;
    }
    q->size++;
}

Customer* dequeue(Queue *q) {
    if (q->head == NULL) return NULL;
    Customer *c = q->head;
    q->head = c->next;
    if (q->head == NULL) q->tail = NULL;
    c->next = NULL;
    q->size--;
    return c;
}

void free_queue(Queue *q) {
    Customer *cur = q->head;
    while (cur) {
        Customer *tmp = cur;
        cur = cur->next;
        free(tmp);
    }
    q->head = q->tail = NULL;
    q->size = 0;
}

/* -------------------------
   Dynamic array for waits
   -------------------------*/

typedef struct DynArray {
    double *data;
    size_t size;
    size_t capacity;
} DynArray;

void dyn_init(DynArray *a) {
    a->size = 0;
    a->capacity = 128;
    a->data = malloc(a->capacity * sizeof(double));
    if (!a->data) { perror("malloc"); exit(1); }
}

void dyn_push(DynArray *a, double v) {
    if (a->size >= a->capacity) {
        a->capacity *= 2;
        a->data = realloc(a->data, a->capacity * sizeof(double));
        if (!a->data) { perror("realloc"); exit(1); }
    }
    a->data[a->size++] = v;
}

void dyn_free(DynArray *a) {
    free(a->data);
    a->data = NULL;
    a->size = a->capacity = 0;
}

/* -------------------------
   Utility math/stat funcs
   -------------------------*/

// comparator for qsort (double)
int cmp_double(const void *pa, const void *pb) {
    double a = *(const double*)pa;
    double b = *(const double*)pb;
    if (a < b) return -1;
    if (a > b) return 1;
    return 0;
}

double compute_mean(DynArray *a) {
    if (a->size == 0) return 0.0;
    double sum = 0;
    for (size_t i=0;i<a->size;i++) sum += a->data[i];
    return sum / (double)a->size;
}

double compute_median(DynArray *a) {
    if (a->size == 0) return 0.0;
    // sort a copy
    double *tmp = malloc(a->size * sizeof(double));
    memcpy(tmp, a->data, a->size * sizeof(double));
    qsort(tmp, a->size, sizeof(double), cmp_double);
    double med;
    if (a->size % 2 == 1) {
        med = tmp[a->size/2];
    } else {
        med = (tmp[(a->size/2)-1] + tmp[a->size/2]) / 2.0;
    }
    free(tmp);
    return med;
}

double compute_stddev(DynArray *a, double mean) {
    if (a->size <= 1) return 0.0;
    double sumsq = 0.0;
    for (size_t i=0;i<a->size;i++) {
        double d = a->data[i] - mean;
        sumsq += d*d;
    }
    return sqrt(sumsq / (double)(a->size)); // population SD; use (n-1) for sample SD
}

int compute_mode_ints(DynArray *a) {
    if (a->size == 0) return 0;
    // mode on integer wait times: we round to nearest integer (or floor)
    // find max wait value to build frequency
    int maxv = 0;
    for (size_t i=0;i<a->size;i++) {
        int v = (int)round(a->data[i]);
        if (v > maxv) maxv = v;
    }
    int *freq = calloc(maxv+1+1, sizeof(int)); // +1 safety
    for (size_t i=0;i<a->size;i++) {
        int v = (int)round(a->data[i]);
        if (v < 0) v = 0;
        freq[v]++;
    }
    int mode = 0;
    int best = -1;
    for (int i=0;i<=maxv;i++) {
        if (freq[i] > best) {
            best = freq[i];
            mode = i;
        }
    }
    free(freq);
    return mode;
}

/* -------------------------
   Random generators
   -------------------------*/

// Return Poisson-distributed integer using Knuth's algorithm
int poisson_rand(double lambda) {
    if (lambda <= 0) return 0;
    double L = exp(-lambda);
    int k = 0;
    double p = 1.0;
    do {
        k++;
        double u = (double)rand() / (RAND_MAX + 1.0);
        p *= u;
    } while (p > L);
    return k - 1;
}

// Return integer service duration, here 2 or 3 minutes with equal probability
int random_service_duration() {
    int r = rand() % 2; // 0 or 1
    return 2 + r; // 2 or 3
}

/* -------------------------
   Main simulation
   -------------------------*/

int main() {
    srand((unsigned int)time(NULL));

    printf("=== Bank Queue (Poisson) Simulator ===\n");
    printf("Simulate an 8-hour bank day (480 minutes).\n");

    double lambda;
    int tellers;
    char cont_after_close_choice;

    printf("Enter average number of customers arriving PER MINUTE (lambda, e.g., 0.2): ");
    if (scanf("%lf", &lambda) != 1) {
        fprintf(stderr, "Invalid input. Exiting.\n");
        return 1;
    }
    printf("Enter number of tellers (e.g., 1): ");
    if (scanf("%d", &tellers) != 1 || tellers < 1) {
        fprintf(stderr, "Invalid number of tellers. Exiting.\n");
        return 1;
    }
    printf("Serve remaining queue after 480 minutes? (y/n) [default y]: ");
    // consume leftover newline
    getchar();
    cont_after_close_choice = getchar();
    int cont_after_close = 1;
    if (cont_after_close_choice == 'n' || cont_after_close_choice == 'N') cont_after_close = 0;

    // State for tellers: each teller has remaining_service_time and pointer to current customer (optional)
    int *teller_busy_remaining = malloc(tellers * sizeof(int));
    Customer **teller_current = malloc(tellers * sizeof(Customer*));
    for (int i=0;i<tellers;i++) { teller_busy_remaining[i] = 0; teller_current[i] = NULL; }

    Queue q;
    init_queue(&q);

    DynArray wait_times;
    dyn_init(&wait_times);

    int CLOCK = 0; // minutes 0..479 for arrivals
    const int END_MINUTE = 480;

    // We'll keep serving until CLOCK reaches 480 and (optionally) queue empty and all tellers idle
    while (1) {
        // Arrivals only for CLOCK < END_MINUTE
        if (CLOCK < END_MINUTE) {
            int arrivals = poisson_rand(lambda);
            for (int a=0;a<arrivals;a++) {
                Customer *c = malloc(sizeof(Customer));
                c->arrival_time = CLOCK;
                c->service_start = -1;
                c->service_duration = 0;
                c->next = NULL;
                enqueue(&q, c);
            }
        }

        // For each teller, if they are idle and queue not empty, start serving next customer
        for (int t = 0; t < tellers; t++) {
            if (teller_busy_remaining[t] <= 0) {
                // if currently serving a customer pointer exists (should not), free it (handled when finished)
                if (teller_current[t] == NULL) {
                    // take next from queue if exists
                    Customer *next = dequeue(&q);
                    if (next != NULL) {
                        next->service_start = CLOCK;
                        next->service_duration = random_service_duration();
                        teller_busy_remaining[t] = next->service_duration;
                        teller_current[t] = next;
                    } else {
                        teller_busy_remaining[t] = 0;
                        teller_current[t] = NULL;
                    }
                } else {
                    // shouldn't happen, but protectively free current if marked done
                    free(teller_current[t]);
                    teller_current[t] = NULL;
                }
            }
        }

        // Advance time by 1 minute: all busy tellers decrement remaining time
        // If a teller finishes this minute (remaining becomes 0), record wait time and free customer
        for (int t = 0; t < tellers; t++) {
            if (teller_busy_remaining[t] > 0) {
                teller_busy_remaining[t]--;
                if (teller_busy_remaining[t] == 0 && teller_current[t] != NULL) {
                    // service finished at minute CLOCK (after decrement)
                    Customer *finished = teller_current[t];
                    int wait = finished->service_start - finished->arrival_time;
                    if (wait < 0) wait = 0;
                    dyn_push(&wait_times, (double)wait);
                    // free finished customer
                    free(finished);
                    teller_current[t] = NULL;
                }
            }
        }

        CLOCK++;

        // termination conditions:
        // - if we've simulated arrival minutes and either:
        //   a) we are NOT serving after close => stop immediately when CLOCK==END_MINUTE
        //   b) we ARE serving after close => stop when CLOCK >= END_MINUTE AND queue empty AND all tellers idle
        if (CLOCK >= END_MINUTE) {
            if (!cont_after_close) {
                // stop right away
                break;
            } else {
                int any_busy = 0;
                for (int t=0;t<tellers;t++) if (teller_busy_remaining[t] > 0) any_busy = 1;
                if (q.size == 0 && !any_busy) break;
                // else continue loop (no more arrivals but serve until done)
            }
        }
        // safety: avoid infinite loop - cap at some large number (e.g., 10000 minutes)
        if (CLOCK > 10000) { fprintf(stderr, "Simulation exceeded safety cap. Breaking.\n"); break; }
    }

    // If some tellers are still serving at exactly break point (e.g., cont_after_close==0 and finishing mid-service),
    // their customers were not counted as finished; depending on brief we ended after 480 minutes. That matches choice.

    // Print basic stats and analysis
    printf("\n--- Simulation finished ---\n");
    printf("Sim minutes simulated (clock): %d\n", CLOCK);
    printf("Customers served (recorded wait times): %zu\n", wait_times.size);
    printf("Customers left in queue: %d\n", q.size);

    if (wait_times.size == 0) {
        printf("No customers finished service; no stats to show.\n");
    } else {
        double mean = compute_mean(&wait_times);
        double median = compute_median(&wait_times);
        double stddev = compute_stddev(&wait_times, mean);
        double maxwait = 0.0;
        for (size_t i=0;i<wait_times.size;i++) if (wait_times.data[i] > maxwait) maxwait = wait_times.data[i];
        int mode = compute_mode_ints(&wait_times);

        printf("\n-- Wait Time Analysis (in minutes) --\n");
        printf("Mean wait time: %.3f\n", mean);
        printf("Median wait time: %.3f\n", median);
        printf("Mode (rounded) wait time: %d\n", mode);
        printf("Standard deviation: %.3f\n", stddev);
        printf("Longest wait time: %.3f\n", maxwait);

        // Optional: show distribution top few values
        printf("\nSample of wait times (first 20): ");
        for (size_t i=0;i < wait_times.size && i < 20; i++) {
            printf("%.0f ", wait_times.data[i]);
        }
        printf("\n");
    }

    // cleanup
    dyn_free(&wait_times);
    for (int t=0;t<tellers;t++) if (teller_current[t] != NULL) free(teller_current[t]);
    free(teller_current);
    free(teller_busy_remaining);
    free_queue(&q);

    printf("\nSave results and add README.md to your repo as required.\n");
    printf("Done.\n");
    return 0;
}
