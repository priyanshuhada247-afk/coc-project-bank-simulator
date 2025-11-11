# coc-project-bank-simulator

**Bank Queue (Poisson) Simulator** — Clash of Coders Capstone Project

## Project Summary
This C program simulates a bank queue over an 8-hour working day (480 minutes). Customer arrivals are modeled using a Poisson distribution (λ = average customers per minute). Customers join a dynamic queue implemented via a linked list (structs + malloc). One or more tellers serve customers; service time is randomized (2–3 minutes). The simulator records each customer's wait time and produces statistical analysis (mean, median, mode, standard deviation, longest wait).

## Concepts Used
**C Concepts**
- `struct` for Customer nodes
- Linked List queue implemented using pointers, `malloc`, and `free`
- Functions for modularity
- `for` loops, `if/else`, dynamic arrays using `realloc`

**Math Concepts**
- Poisson distribution (Knuth's algorithm)
- Central tendencies: Mean, Median, Mode
- Standard deviation

## Files
- `bank_simulator.c` — main C source file

## How to Compile
```bash
gcc bank_simulator.c -o bank_simulator -lm
