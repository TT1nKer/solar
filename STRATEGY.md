# Strategic Roadmap: From Orbital Engine to Universe Operating System

> Status: Phases 1-10 complete (orbital mechanics core). This document defines the next evolution.

## Current State Assessment

The orbital engine is solid:
- Celestial mechanics (Kepler, N-body, Lambert, CR3BP, Halo orbits)
- Precision ephemeris (DE440, sub-meter)
- Modular physics (13 force models composable)
- Adaptive integration (DOPRI5 + Verlet)
- Mission simulation (Mars Direct, Lunar Gateway, Voyager Grand Tour)

**What's strong**: "how the universe moves" — dynamical correctness.

**What's missing**: "should we fly?", "what if it fails?", "where does supply come from?", "how do multiple agents coordinate?", "how to keep expanding after comms break?"

The limiting factor is no longer numerical accuracy. It's **engineering closure and autonomous closure**.

## Architecture: 8-Layer Universe Operating System

### Layer 1: Celestial Core (DONE)
Ephemeris, time systems, coordinate transforms, constants, reference frames.
- Refactor into fully independent module with unified state interface
- All higher layers read through this — single source of truth

### Layer 2: Dynamics Core (DONE)
N-body propagation, perturbation models, Lambert, flyby, CR3BP, stability.
- Split into `predict()` (propagate) and `design()` (inverse-solve) interfaces
- Future agents call interfaces, never touch internals

### Layer 3: Vehicle Layer (NEXT PRIORITY)
Real spacecraft are not point masses.

**Vehicle state**: dry mass, propellant, battery/SOC, power generation, thermal, attitude control, comm bandwidth, component health, compute resources, spare parts.

**Subsystem models**: propulsion (thrust limits, ignition count, failure rate), power (solar/RTG/battery), thermal (limits, derating), comms (distance-rate relation, occultation, delay), ADCS (pointing cost, constraints), payload (observation/mining/manufacturing efficiency).

This is where the system stops being "point mass trajectories" and becomes "degradable engineering entities."

### Layer 4: Uncertainty Layer (HIGHEST VALUE)
Every result becomes: nominal value + variance + confidence + risk label.

- **Error propagation**: initial state, thrust, mass, timing, model errors
- **Monte Carlo**: "82.4% reachable in 1000 samples; failure from narrow capture window; most sensitive to thrust bias"
- **Sensitivity analysis**: which parameter matters most, where to spend sensor budget, which phases need replanning

The most important capability for AI in space: knowing **where it's most likely wrong**.

### Layer 5: Infrastructure & Logistics
Stop thinking "one flight." Start thinking "network."

**Nodes**: planetary orbits, Lagrange points, asteroid stations, relay stations, mining sites, manufacturing sites, depots.

**Edges**: Lambert transfers, low-energy transfers, periodic launch windows, supply lines, comm visibility.

**Edge weights**: delta-v, time, risk, fuel, comm delay, success probability, transport capacity.

This answers: which node to occupy first, which supply line is most stable, where to build long-term, which edges are cargo vs crew.

### Layer 6: Resource & Industry
Abstract resource flow (before detailed chemistry):
- asteroid → ore
- ore → metal + waste (refinery)
- metal + energy → spare_parts (factory)
- water + energy → propellant (depot)

Once there's resource flow, the system goes from "Earth sends a ship" to "system can self-sustain in the environment." This is where expansion actually begins.

### Layer 7: Agent Layer
Not just "user writes plan, system executes."

**Agent interface**: `observe()`, `estimate()`, `plan()`, `act()`, `replan()`, `negotiate()`.

**Agent types**: spacecraft agent, relay station agent, mining agent, manufacturing agent, mission control agent, civilization strategy agent.

**Implementation order**: rule-based → planner (A*/MPC/tree search) → learning (RL/imitation) → LLM integration last.

### Layer 8: Civilization Strategy
Optimize civilization-level metrics, not single trajectories:
- Maximum sustainable years
- Maximum node survival redundancy  
- Minimum Earth dependency
- Maximum expansion efficiency per unit resource
- Maximum knowledge replication rate
- Maximum survival probability after communication loss

Strategic modes: exploration-first, resource-first, redundancy-first, expansion-first, survival-first.

## Development Order

### Next-1: Vehicle + Uncertainty
- Spacecraft subsystems (propulsion, power, comms, health)
- Monte Carlo runner
- Sensitivity analysis
- Fault injection

**Transforms project from**: orbital solver → **mission simulator**

### Next-2: Network + Logistics  
- Node graph with solar system topology
- Edge weight generator (delta-v, time, risk)
- Multi-objective shortest path
- Depot and cargo flow

**Transforms project from**: single mission chain → **solar system transportation network**

### Next-3: Resource + Manufacturing
- Abstract resource model
- Mining/refining/manufacturing pipeline
- Inventory and transport
- Node production capacity

**Transforms project from**: transportation system → **preliminary space economy**

### Next-4: Autonomous Agents
- Rule-based agents
- Planner agents  
- Learning agents (later)

**Transforms project from**: human plans, machine executes → **machine plans, replans, cooperates**

### Next-5: Civilization Sandbox
Strategy layer. This is where "the universe system" begins.

## Anti-Hollowing Protocol

### Rules for AI-driven development:
1. Task unit = **capability closure**, not module
2. Every delivery requires: code + test + sample run + unimplemented list
3. No structure without behavior evidence
4. Vertical slices before horizontal infrastructure
5. Tests first, implementation second
6. Adversarial self-audit after each delivery
7. 20% architecture, 80% runnable slices
8. No file created unless called in current execution chain
9. No placeholder returns disguised as completion
10. README is the last thing to trust — sample run is the first

### Validation template for each delivery:
```
ACTUALLY COMPLETED:
- [specific runnable capability with evidence]

NOT COMPLETED:
- [explicit list, not hidden behind "extensible"]

MOST LIKELY BUGS:
- [honest assessment]

FASTEST WAY TO FALSIFY:
- [specific input → expected output, if wrong = broken]
```

## Data Structure Evolution

### Separate concerns:
- `CelestialBody`: mass, radius, rotation, ephemeris, gravity field, atmosphere, resources
- `Vehicle`: state vector, subsystems, resource inventory, mission queue, health, agent
- `InfrastructureNode`: orbit/position, storage, production, comms, redundancy

### Mission evolution:
From event lists → goal + constraints + policies + contingencies:
```
Goal: deliver 2000 kg propellant to Mars transfer depot
Constraints: failure probability < 5%, arrival before 2035-06-01
Policy: prefer reusable infrastructure
Contingencies: if window missed, reroute via staging node
```

### All modules output explainable results:
- Why this path
- Cost breakdown
- Risk source
- Most sensitive parameter
- Recommended replan on failure
