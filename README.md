# S-CFLAP_Code
Java research code for a sequential competitive facility location and attractiveness problem. The project models leader-follower competition, selecting facility locations and attractiveness levels to optimize market share.

## Project Structure

- `src/BC/`: Continuous-attractiveness models, including outer approximation (OA), second-order cone programming (SOCP), and leader-follower solution algorithms.
- `src/BC_Discrete/`: Models and algorithms with discrete attractiveness levels.
- `src/Heuristic/`: PEH and particle swarm optimization (PSO) implementations.
- `src/Common/`: Shared instance, facility, customer, solution, and cut classes.
- `RandomData/`: Experiment input location and instance-name list.
- `RandomResults/`: Experiment output location.

## Requirements

- Java 8 (the version configured in the Eclipse project).
- IBM ILOG CPLEX with its Java API and native libraries. The existing configuration references CPLEX Studio 12.8.
- Eclipse IDE, or an equivalent Java development environment.
