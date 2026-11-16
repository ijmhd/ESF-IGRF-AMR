# Runtime function reference

Every runtime C definition is listed below; see the [complete annotated-source atlas](../methodology/source-atlas.pdf).

## `src/amr3.c`

AMR driver and conservative O+ fluxes. n is normalized by 1e12 m^-3; m=n*u_parallel. Private levels -1/0/1 represent physical L0/L1/L2 in a three-level run.

* [`wrap`](../src/amr3.c#L18): Return the nonnegative periodic index i modulo n.
* [`cover`](../src/amr3.c#L26): Map a canonical transverse cell to its backing leaf column; apply the configured periodic/clamped indexing.
* [`owned`](../src/amr3.c#L39): Test whether a backing column belongs to this rank, so shared faces update only local conserved states.
* [`ph`](../src/amr3.c#L45): Read a transverse potential node using this solver's boundary/periodic indexing convention.
* [`derivative`](../src/amr3.c#L52): Compute the two transverse potential derivatives at a geometric face and retain both metric cross contributions.
* [`minmod`](../src/amr3.c#L78): Return the sign-consistent slope with smallest magnitude; return zero at an extremum.
* [`length_at`](../src/amr3.c#L84): Read the physical along-field segment length at the requested private refinement level.
* [`cell_value`](../src/amr3.c#L103): Sample density or velocity with the caller's parallel boundary convention; used by the along-field reconstruction.
* [`parallel_state_at`](../src/amr3.c#L112): Construct a bounded left/right parallel face state at an explicitly requested refinement level.
* [`parallel_state`](../src/amr3.c#L122): Resolve the actual local density level, then reconstruct its parallel face state.
* [`length`](../src/amr3.c#L129): Resolve a backing cell's active density level before reading its physical parallel length.
* [`log_pressure`](../src/amr3.c#L137): Compute the nonuniform parallel logarithmic pressure gradient, using common Ne in the multi-ion closure and exact-vacuum secants.
* [`bernoulli`](../src/amr3.c#L167): Evaluate the Bernoulli fitting function z/(exp(z)-1) with small-argument series and large-argument limits.
* [`parallel_cache`](../src/amr3.c#L182): Prepare shared parallel flux/rate tables for a column-uniform density layout; segmented layouts use a different path.
* [`density_cover`](../src/amr3.c#L246): Locate the actual density owner at a canonical cell and along-field index; owner can differ from electric-field backing.
* [`append_face`](../src/amr3.c#L260): Append or merge a canonical face interval whose two active density owners and levels are unchanged.
* [`a3_face_topology`](../src/amr3.c#L299): Build conservative face blocks for the current segmented density mesh; identical interface fluxes are shared by both adjacent volumes.
* [`metric_at`](../src/amr3.c#L332): Read or conservatively aggregate geometry for a requested physical density level, including the macro level.
* [`a3_parallel_geometry`](../src/amr3.c#L348): Cache static segment lengths and parallel relaxation factors for every supported level.
* [`segmented_parallel`](../src/amr3.c#L389): Evaluate parallel particle/momentum flux for one actual segmented interface; retain collisional relaxation and electron-pressure fitting.
* [`marker_fraction`](../src/amr3.c#L453): Return a bounded material seed fraction from the advected marker and density; zero density gives zero fraction.
* [`a3_transport`](../src/amr3.c#L464): Compute O+ finite-volume density/momentum RHS, shared face fluxes, gravity/pressure sources and global stability/travel rates. Collective; sampling may populate lazy caches, so AMR face loops remain on the main thread.
* [`rhs`](../src/amr3.c#L640): Execute the complete current-stage chain: prepare/halo, integrated potential solve, O+ transport, then molecular transport. The potential always belongs to this stage's density and chemistry state.
* [`node_id`](../src/amr3.c#L652): Flatten node coordinates with the along-field coordinate contiguous.
* [`make_box`](../src/amr3.c#L659): Build Cartesian synthetic geometry for transport/potential verification rather than tracing a physical magnetic field.
* [`gather_geometry`](../src/amr3.c#L700): Collect disjoint computed geometry row ranges into replicated canonical geometry, in MPI-count-safe chunks.
* [`parallel_geometry`](../src/amr3.c#L717): Distribute tube tracing and cell/face metric construction among ranks, then assemble the replicated geometry and diagnostics.
* [`geometry`](../src/amr3.c#L772): Load and validate a geometry cache, or construct geometry and optionally write a new cache; initialize conductivity responses afterward.
* [`output`](../src/amr3.c#L833): Open a named output in the current run directory; failures are collective-fatal through a3_fail.
* [`metadata`](../src/amr3.c#L845): Write effective geometry, numerical controls, units and binary layout metadata for this AMR run.
* [`snapshot`](../src/amr3.c#L920): Collect owned states, reconstruct canonical display density, and write time/mesh-bound outputs plus the native checkpoint. Display reconstruction is not an evolved finest-grid solution.
* [`remap_tests`](../src/amr3.c#L1105): Exercise constant/background/smooth conservative transfers and validate particle/momentum invariants in the dedicated remap case.
* [`import_static`](../src/amr3.c#L1198): Import a uniform-grid restart by preserving original coarse-cell particle and signed momentum integrals on the canonical AMR geometry.
* [`workspace`](../src/amr3.c#L1285): Replace scratch arrays after a mesh rebuild so their sizes match the new owned-plus-halo backing layout.
* [`amr3_run`](../src/amr3.c#L1296): Drive initialization/restart, regridding, half chemistry, two-stage O+/molecular transport, acceptance checks, diagnostics and checkpointing. All ranks follow the same clocks; no level subcycling is used.

## `src/amr3_adapt.c`

Collective mesh rebuild, ownership migration and conservative transfer. Geometry remains fixed; density and electric-potential approximation spaces may change separately.

* [`migrate`](../src/amr3_adapt.c#L16): Transfer O+ density/momentum/optional marker intersections between old and new rank-owned slabs before conservative remapping.
* [`a3_adapt`](../src/amr3_adapt.c#L47): Detect actual mask/field-space changes, repartition molecules and O+, build/remap the new mesh, verify global invariants and resolve potential. A no-change request returns without rebuilding.

## `src/amr3_cache.c`

Native-ABI geometry cache. This is a performance cache, not a portable state checkpoint. Source/input identities, ABI sizes and payload checksum must all match.

* [`hash`](../src/amr3_cache.c#L20): Update the FNV-1a checksum over an exact byte sequence.
* [`file_hash`](../src/amr3_cache.c#L31): Hash an input file for cache identity; an empty optional path has identity zero.
* [`a3_geometry_identity`](../src/amr3_cache.c#L51): Accept the exact compiled geometry identity or one explicitly enumerated, version-guarded historical equivalence; reject unknown changes.
* [`a3_geometry_cache`](../src/amr3_cache.c#L88): Read/write native canonical geometry with complete parameter/input/source/ABI checks, atomic rename and payload checksum.

## `src/amr3_checkpoint2.c`

Structure-aware ESFAMR02/03/04/05 checkpoint formats. Preserve density/field masks, hysteresis clocks, optional seed tracer and conservation history.

* [`bytes`](../src/amr3_checkpoint2.c#L25): Read or write exactly the requested bytes while updating the checkpoint checksum; short transfers are fatal.
* [`integer`](../src/amr3_checkpoint2.c#L38): Encode/decode an unsigned 64-bit integer in little-endian order and update the stream checksum.
* [`reals`](../src/amr3_checkpoint2.c#L52): Transfer floating data in portable little-endian order, retaining the checksum over the actual encoded bytes.
* [`finish`](../src/amr3_checkpoint2.c#L68): Write or verify the final checksum, reject trailing data on reads, and close the checkpoint stream.
* [`parameters`](../src/amr3_checkpoint2.c#L80): Choose the exact restart signature length required by this case's driver/neutral extensions.
* [`signature`](../src/amr3_checkpoint2.c#L87): Build the ordered physical/numerical restart signature; changing a signed parameter intentionally rejects incompatible restart.
* [`a3_checkpoint2_write`](../src/amr3_checkpoint2.c#L127): Atomically write structure-aware masks/history, optional tracer, actual conserved states, potential and checksummed physical/numerical signatures.
* [`a3_checkpoint2_read`](../src/amr3_checkpoint2.c#L209): Validate and reconstruct a structure-aware checkpoint, distribute states under the current MPI decomposition, and restore mesh/clocks/history.

## `src/amr3_macro.c`

Third physical level. A macro cell covers four internal roots. Only the primary backing slot evolves when its macro mask is coarse; volume-weighted virtual children serve reconstruction.

* [`a3_nx`](../src/amr3_macro.c#L18): Return the transverse x dimension for a private level; level -1 is the extra physical coarse level.
* [`a3_ny`](../src/amr3_macro.c#L25): Return the transverse y dimension for a private level; level -1 is the extra physical coarse level.
* [`a3_scale`](../src/amr3_macro.c#L31): Return the number of canonical finest cells across one requested-level cell: 4, 2 or 1.
* [`a3_macro_index`](../src/amr3_macro.c#L37): Map an internal root index to its containing two-by-two macro cell.
* [`a3_macro_primary`](../src/amr3_macro.c#L43): Return the unique backing column used to evolve an unrefined macro cell.
* [`a3_rebase`](../src/amr3_macro.c#L50): Repair interior pointers after copying an A3 struct; they must point into the destination's own embedded storage.
* [`coord`](../src/amr3_macro.c#L59): Wrap or clamp macro coordinates according to physical versus synthetic boundary policy.
* [`mm`](../src/amr3_macro.c#L69): Apply the two-argument minmod slope limiter without changing a conserved cell average.
* [`present`](../src/amr3_macro.c#L75): Test whether this rank holds the primary backing state for a macro cell.
* [`around`](../src/amr3_macro.c#L82): Read a neighboring macro representation, falling back to the central cell if the neighbor is outside the local halo.
* [`a3_macro_init`](../src/amr3_macro.c#L95): Allocate macro masks/history/state tables and copy their history from the old mesh; reuse immutable volume/background values on rebuild.
* [`a3_macro_metrics`](../src/amr3_macro.c#L159): Reconstruct unsegmented raw metrics from fixed canonical cells before applying any density mask. Repeated invocation must not aggregate already aggregated volumes.
* [`a3_macro_prepare`](../src/amr3_macro.c#L190): Restrict evolved owner states to macro averages, then construct bounded conservative virtual internal-root states for coarse macro cells.
* [`a3_macro_sample`](../src/amr3_macro.c#L336): Read a prepared macro density/momentum/velocity/marker sample and verify local halo availability.
* [`a3_macro_tags`](../src/amr3_macro.c#L348): Choose coarse-level refinement from residual gradient/curvature and parent information loss; enforce nesting, buffering and physical-time quiet hysteresis.
* [`a3_macro_output`](../src/amr3_macro.c#L454): Write physical L0/L1/L2 cell labels and actual activity counts for visualization and validation.
* [`node`](../src/amr3_macro.c#L487): Flatten a transverse potential node index on the canonical nodal grid.
* [`edge`](../src/amr3_macro.c#L495): Determine a macro-edge constraint in the transverse potential approximation space.
* [`add_map`](../src/amr3_macro.c#L502): Accumulate constrained basis weights for one potential map without losing contributions shared through a parent.
* [`resolve`](../src/amr3_macro.c#L518): Resolve nested potential-node constraints to the final independent DOFs.
* [`a3_macro_maps`](../src/amr3_macro.c#L549): Build nested bilinear potential constraints for the extra coarse level, keeping boundary and hanging-node consistency.

## `src/amr3_mesh.c`

Mesh, owned/ghost storage and conservative reconstruction. Backing columns are not all active density DOFs. Virtual samples may lazily populate caches; this module is not generally thread-safe.

* [`a3_fail`](../src/amr3_mesh.c#L18): Report rank and invariant violation, abort the communicator, and never continue with an invalid model state.
* [`a3_alloc`](../src/amr3_mesh.c#L27): Allocate zero-initialized storage, including a nonnull minimum allocation for empty logical arrays; allocation failure is fatal.
* [`wrap`](../src/amr3_mesh.c#L36): Return the nonnegative periodic index i modulo n.
* [`index_axis`](../src/amr3_mesh.c#L43): Wrap synthetic/periodic indices and otherwise clamp physical boundary indices.
* [`minmod`](../src/amr3_mesh.c#L51): Return the sign-consistent slope with smallest magnitude; return zero at an extremum.
* [`gid`](../src/amr3_mesh.c#L57): Flatten a density-level cell index; k is contiguous and nx depends on the private level.
* [`moment_weights`](../src/amr3_mesh.c#L65): Integrate canonical physical-volume moments, restrict them to coarser levels, and build four-cell face weights for PDM4.
* [`geometry_cache`](../src/amr3_mesh.c#L168): Share reference-counted immutable geometry/response/PDM tables across regrids and allocate the new metric table. Three-level raw metrics are built once by a3_macro_metrics.
* [`a3_partition_plan`](../src/amr3_mesh.c#L278): Balance the existing density-work estimate plus canonical-field weight using whole transverse root rows. Three-level cuts align to pairs; this one-dimensional policy has a finite rank limit.
* [`a3_mesh`](../src/amr3_mesh.c#L320): Construct a mesh from the requested fixed/uniform/dynamic mode without an explicit supplied mask.
* [`a3_mesh_mask`](../src/amr3_mesh.c#L327): Construct a mesh with explicit density masks and optional old-mesh geometry reuse / ownership cuts.
* [`a3_mesh_mask_work`](../src/amr3_mesh.c#L335): Construct backing columns, rank ownership, halos, reconstruction workspaces, geometry metrics and face topology from explicit masks. Collective layout decisions must agree on all ranks.
* [`a3_segment_metrics`](../src/amr3_mesh.c#L550): Reset and aggregate physical metrics for the actual segmented density owners; inactive backing slots receive zero volume.
* [`exchange`](../src/amr3_mesh.c#L616): Exchange owned O+ density/momentum and optional marker halos using nonblocking peer transfers, then unpack only after completion.
* [`a3_prepare`](../src/amr3_mesh.c#L655): Bind the current stage, exchange state, reset lazy caches, restrict actual conserved states, construct macro representations and refresh molecular/electron responses. Must finish before potential or flux sampling.
* [`parent`](../src/amr3_mesh.c#L715): Read a parent density/momentum sample after applying boundary indices; reject a stencil outside the exchanged root-row halo.
* [`child_proposal`](../src/amr3_mesh.c#L728): Evaluate a conservative background-shaped density child from cached parent amplitude and residual slopes.
* [`slopes`](../src/amr3_mesh.c#L742): Build bounded background-aware density reconstruction and lazy child-density caches for one parent and along-field index.
* [`velocity_slopes`](../src/amr3_mesh.c#L800): Prolong parallel velocity with particle-weighted centering so child momentum sums to the parent momentum.
* [`a3_proposed_child`](../src/amr3_mesh.c#L846): Return a prospective child density/momentum for transfer/coarsening assessment while preserving the existing reconstruction semantics.
* [`a3_sample`](../src/amr3_mesh.c#L863): Sample actual or virtual states at a requested level. It may update lazy density/velocity child caches; concurrent calls require explicit cache ownership.
* [`sample_pair`](../src/amr3_mesh.c#L902): Collect a one-dimensional face stencil at a common requested level, including virtual samples across coarse/fine contacts.
* [`a3_face_pair`](../src/amr3_mesh.c#L938): Return particle-density and velocity face states using physical-volume PDM4 or the configured MUSCL fallback.
* [`initial_fine`](../src/amr3_mesh.c#L977): Evaluate the idealized Chapman/seed fixture on canonical physical geometry; used for fresh initialization and verification.
* [`a3_initialize_marker`](../src/amr3_mesh.c#L1027): Initialize an advected seed-protection material quantity consistently with density; this tracer influences tagging, not conductivity.
* [`a3_initialize`](../src/amr3_mesh.c#L1067): Restrict the canonical idealized state to actual active volumes; leave nonowner backing states inactive.
* [`a3_mass`](../src/amr3_mesh.c#L1104): Integrate a conserved field with actual owner volumes and sum across MPI ranks.
* [`a3_valid`](../src/amr3_mesh.c#L1116): Collectively reject nonfinite/negative density, invalid vacuum momentum and out-of-range marker states; never impose a density floor.
* [`a3_remap`](../src/amr3_mesh.c#L1132): Conservatively restrict/prolong density, parallel momentum and optional marker between old/new active meshes, checking global transfer invariants.
* [`a3_free`](../src/amr3_mesh.c#L1182): Release per-mesh arrays and PETSc objects; free immutable shared caches only when their reference count reaches zero.

## `src/amr3_migration.c`

Molecular old/new slab ownership intersections. Separate receive buffers preserve every old owned value until O+ remapping and all sends finish.

* [`a3_night_repartition`](../src/amr3_migration.c#L20): Move only molecular cells whose owner changes using old/new slab intersections and separate receive storage. Identical intervals are a byte-preserving no-op; old owned values survive until O+ remap completes.

## `src/amr3_molecular.c`

NO+ and O2+ transport on canonical finest cells, with separate parallel momenta and common electron pressure. Canonical geometry is fixed even when O+ adapts.

* [`clamp`](../src/amr3_molecular.c#L20): Clamp a cell index to the nearest valid boundary cell.
* [`id`](../src/amr3_molecular.c#L26): Flatten a canonical molecular cell with periodic-x or clamped physical indices.
* [`nid`](../src/amr3_molecular.c#L37): Flatten node coordinates with the along-field coordinate contiguous.
* [`center`](../src/amr3_molecular.c#L43): Average the four geometric vertices of a canonical face in ECEF metres.
* [`a3_molecular_geometry`](../src/amr3_molecular.c#L59): Build owned canonical face mobility/gravity coefficients, and rebuild only when the owned slab changes.
* [`a3_molecular_workspace`](../src/amr3_molecular.c#L137): Resize owned molecular stage, backup, rate, pressure-log and same-attempt ETD work arrays after an ownership change.
* [`a3_night_electrons`](../src/amr3_molecular.c#L175): Reconstruct current canonical O+ plus molecules to form quasineutral Ne and invalidate every averaged-Ne cache; lazy O+ sampling keeps this on the main thread.
* [`a3_night_electron`](../src/amr3_molecular.c#L197): Return the original volume-weighted current-stage Ne average, memoized only for owned cells and invalidated whenever electrons are rebuilt. Main-thread-only lazy cache; no density floor or altered summation.
* [`length`](../src/amr3_molecular.c#L224): Read the cached canonical physical parallel segment length.
* [`cached_jump`](../src/amr3_molecular.c#L231): Use current-stage cached logarithms for positive states and the original bounded secant at exact vacuum.
* [`slope`](../src/amr3_molecular.c#L239): Compute a nonuniform centered/one-sided parallel logarithmic gradient from the current-stage pressure cache.
* [`minmod`](../src/amr3_molecular.c#L256): Return the sign-consistent slope with smallest magnitude; return zero at an extremum.
* [`state`](../src/amr3_molecular.c#L263): Reconstruct a bounded parallel molecular density or velocity face state using the three-point slope limiter.
* [`a3_night_transport`](../src/amr3_molecular.c#L278): Evaluate both species' conservative canonical fluxes, parallel pressure/gravity sources and stability rates. Threaded sweeps assign disjoint full grid lines; MPI reductions remain on the main thread.
* [`valid`](../src/amr3_molecular.c#L441): Collectively verify finite/nonnegative molecular densities and consistent exact-vacuum momenta.
* [`a3_night_transport_begin`](../src/amr3_molecular.c#L457): Save the accepted pre-transport molecular state and first-stage boundary rate for the two-stage update.
* [`a3_night_predict`](../src/amr3_molecular.c#L472): Apply the molecular first-stage density/ETD momentum predictor and cache its fixed-nu coefficients for this exact attempted dt.
* [`a3_night_correct`](../src/amr3_molecular.c#L497): Complete the two-stage density and ETD2 momentum update using same-attempt coefficients, then check positivity and stage boundary balance.
* [`a3_night_accept`](../src/amr3_molecular.c#L520): Commit the accepted molecular boundary-particle increment after the coupled step passes all checks.

## `src/amr3_night.c`

Experimental isothermal, night-only multi-ion closure. Quasineutral Ne equals O+ + NO+ + O2+. Chemistry uses SI densities internally; stored ion densities are normalized by 1e12 m^-3.

* [`a3_night_loss`](../src/amr3_night.c#L22): Return cumulative O+ chemical conversion loss for the global O+ conservation budget.
* [`gravity_drift`](../src/amr3_night.c#L29): Solve perpendicular collisional gravity mobility with signed gyrofrequency and remove the parallel force component.
* [`halo`](../src/amr3_night.c#L38): Exchange molecular density/momentum ghost slabs; receives and sends remain valid and disjoint until MPI_Waitall.
* [`a3_night_gather`](../src/amr3_night.c#L71): Collect the full canonical molecular state for saved outputs/checkpoints only; normal RHS/regrid paths use halos and owned migration.
* [`a3_night_refresh`](../src/amr3_night.c#L91): Refresh current molecular Pedersen/Hall and gravity-current responses in owned plus ghost cells after halo exchange.
* [`bytes`](../src/amr3_night.c#L137): Read or write exactly the requested bytes while updating the checkpoint checksum; short transfers are fatal.
* [`integer`](../src/amr3_night.c#L149): Encode/decode an unsigned 64-bit integer in little-endian order and update the stream checksum.
* [`reals`](../src/amr3_night.c#L163): Transfer floating data in portable little-endian order, retaining the checksum over the actual encoded bytes.
* [`finish`](../src/amr3_night.c#L179): Write or verify the final checksum, reject trailing data on reads, and close the checkpoint stream.
* [`read_restart`](../src/amr3_night.c#L192): Validate and read the molecular checkpoint, numerical-source identity and binding to the native checkpoint checksum.
* [`a3_night_create`](../src/amr3_night.c#L250): Load frozen neutral/rate inputs, initialize NO+/O2+ reservoirs and transport coefficients, or restore a paired molecular checkpoint; enforce the night-only applicability guard.
* [`a3_night_begin`](../src/amr3_night.c#L378): Back up molecular state and reaction budgets before an attempted coupled step.
* [`a3_night_restore`](../src/amr3_night.c#L395): Restore the pre-attempt molecular state and reaction budgets when a CFL or positivity check requires another attempt.
* [`macro_chemistry`](../src/amr3_night.c#L410): Apply local reactions to disjoint canonical children of actual O+ volumes, then volume-average O+ back to its owner. Thread-local budgets are reduced; fatal handling stays outside workers.
* [`a3_night_step`](../src/amr3_night.c#L469): Execute one chemistry split operator, check its charge budget, reduce global conversion/recombination totals, and refresh conductivity responses.
* [`a3_night_snapshot`](../src/amr3_night.c#L546): Write a molecular checkpoint bound to the just-written native checkpoint, molecular frames, metadata and total-ion reaction/boundary budgets.
* [`a3_night_free`](../src/amr3_night.c#L644): Release all molecular full-grid and owned-work arrays, including stage-local ETD/log caches.

## `src/amr3_potential.c`

Field-line-integrated electrostatic operator on constrained bilinear transverse elements. Both Hall cross terms are retained; phi is constant along the traced field coordinate.

* [`a3_driver_active`](../src/amr3_potential.c#L32): Report whether a prescribed constant or time-dependent external eastward field is active.
* [`a3_driver_field`](../src/amr3_potential.c#L38): Evaluate the external eastward field at the actual SSPRK stage time.
* [`a3_driver_gradient`](../src/amr3_potential.c#L48): Convert eastward field to the model's transverse potential-label gradient using the recorded reference length.
* [`vertex`](../src/amr3_potential.c#L55): Flatten a canonical transverse potential vertex.
* [`boundary`](../src/amr3_potential.c#L64): Classify a potential node or apply its prescribed boundary role; inspect the explicit regional/periodic branches below.
* [`map_add`](../src/amr3_potential.c#L70): Accumulate one constrained-node basis contribution, merging duplicate independent DOFs.
* [`a3_potential_create`](../src/amr3_potential.c#L85): Create constrained field maps, distributed sparse PETSc objects, FGMRES/GAMG configuration and optional precomputed COO basis tables.
* [`coo_create`](../src/amr3_potential.c#L175): Precompute constrained bilinear element bases and PETSc COO sparsity/index maps; material coefficients remain current-stage data.
* [`integrate_tensor`](../src/amr3_potential.c#L289): Integrate both Pedersen/Hall tensor cross terms and gravity-current forcing along the full field line using canonical physical volumes.
* [`a3_potential_seed`](../src/amr3_potential.c#L318): Seed independent potential DOFs from a prior canonical potential for a warm start after regridding.
* [`element`](../src/amr3_potential.c#L333): Assemble one bilinear element with the full nonsymmetric tensor and prescribed potential contributions using the reference element path.
* [`regional_sides`](../src/amr3_potential.c#L391): Impose the regional polarization-potential normal-derivative relation, retaining Hall and nonorthogonal metric terms.
* [`a3_potential`](../src/amr3_potential.c#L431): Assemble the current-stage integrated operator, solve with controlled preconditioner reuse, explicitly verify the true residual, and scatter canonical potential. Collective; no worker-thread PETSc calls.
* [`a3_potential_free`](../src/amr3_potential.c#L514): Destroy constrained maps, COO workspaces and every PETSc object owned by this mesh.

## `src/amr3_restart.c`

Legacy ESFAMR01 portable checkpoint and dispatch to structure-aware versions. Integer and floating data are encoded little-endian and bound to input signatures.

* [`bytes`](../src/amr3_restart.c#L28): Read or write exactly the requested bytes while updating the checkpoint checksum; short transfers are fatal.
* [`integer`](../src/amr3_restart.c#L41): Encode/decode an unsigned 64-bit integer in little-endian order and update the stream checksum.
* [`reals`](../src/amr3_restart.c#L55): Transfer floating data in portable little-endian order, retaining the checksum over the actual encoded bytes.
* [`signature`](../src/amr3_restart.c#L71): Build the ordered physical/numerical restart signature; changing a signed parameter intentionally rejects incompatible restart.
* [`finish`](../src/amr3_restart.c#L86): Write or verify the final checksum, reject trailing data on reads, and close the checkpoint stream.
* [`a3_checkpoint_write`](../src/amr3_restart.c#L99): Write ESFAMR01 for legacy masks or dispatch to the structure-aware writer; retain little-endian encoding and atomic replacement.
* [`a3_checkpoint_read`](../src/amr3_restart.c#L170): Recognize the checkpoint magic, validate legacy signatures or dispatch to the structure reader, then scatter conserved states to current rank ownership.

## `src/amr3_structure.c`

Solution-based refinement detector, including represented-background residuals, velocity curvature, seed protection and prospective-parent information loss. Detector thresholds never modify plasma density or conductivity.

* [`a3_owner`](../src/amr3_structure.c#L31): Return the unique evolved density backing owner for a column and k; account for both segmented L1/L2 and macro L0 ownership.
* [`a3_active`](../src/amr3_structure.c#L40): Test whether this backing slot is the actual evolved density owner at k.
* [`a3_density_leaf`](../src/amr3_structure.c#L47): Return the physical density cell's level and labels, which can differ from the electric-field/backing column level.
* [`a3_density_mask`](../src/amr3_structure.c#L65): Select the segmented density mask or legacy column mask as required by the active case.
* [`index1`](../src/amr3_structure.c#L71): Apply boundary indexing on the detector's root/fine transverse stencil.
* [`drift`](../src/amr3_structure.c#L83): Evaluate or read a cached current E-cross-B drift averaged over canonical children at the requested detector level.
* [`spacing`](../src/amr3_structure.c#L114): Estimate physical transverse cell spacing from contravariant metric magnitude, solely for detector wavelength/velocity normalization.
* [`indicator`](../src/amr3_structure.c#L124): Compute normalized 3x3 density gradient/Hessian and drift-curvature indicators; optionally remove the locally fitted represented-background component.
* [`integrated`](../src/amr3_structure.c#L183): Compute/cache the root-cell field-line-integrated conductivity tensor used to detect elliptic-coefficient curvature.
* [`coefficient_curve`](../src/amr3_structure.c#L203): Measure normalized transverse second differences of the complete integrated tensor, including independent Hall cross terms.
* [`mode_needs_fine`](../src/amr3_structure.c#L223): Check declared seed/envelope sampling requirements and count locations that remain unresolved even on the finest available mesh.
* [`a3_structure_tags`](../src/amr3_structure.c#L240): Collectively form structure-based refinement requests, prospective-parent loss tests, seed guard, slab promotion, buffers and hysteresis. No unperturbed control solution is evolved by this detector.
* [`out`](../src/amr3_structure.c#L426): Open one detector diagnostic output or abort on failure.
* [`a3_structure_output`](../src/amr3_structure.c#L438): Write actual density masks, tag reason bits/values, marker diagnostics and frame-to-mesh bindings for plotting.

## `src/amr3_tags.c`

Legacy percentile detector and dispatch to the structure detector. The published teaching cases explicitly select structure; legacy behavior remains available for comparison.

* [`order`](../src/amr3_tags.c#L17): Compare two finite density values for the legacy longitude-percentile detector.
* [`wrap`](../src/amr3_tags.c#L24): Return the nonnegative periodic index i modulo n.
* [`a3_tags`](../src/amr3_tags.c#L32): Count a tag check and dispatch either the structure criterion or legacy percentile/depletion criterion; both only request mesh changes.

## `src/electric_driver.c`

Validated piecewise-linear imposed eastward electric field. Times are elapsed model seconds; field units in the CSV are V/m.

* [`electric_driver_load`](../src/electric_driver.c#L14): Read/hash a strict time_s,east_V_m CSV, require finite increasing times from zero and coverage through the requested absolute end time.
* [`electric_driver_value`](../src/electric_driver.c#L55): Linearly interpolate the imposed field at a model stage time; return NaN outside declared coverage.
* [`electric_driver_next`](../src/electric_driver.c#L77): Return the next forcing knot after time, or infinity if absent, so accepted steps land on forcing changes.

## `src/field_cli.c`

Small IGRF command-line diagnostic; angles are degrees, positions metres and magnetic output nT.

* [`main`](../src/field_cli.c#L13): Parse this executable's interface, validate inputs, run its diagnostic/model entry point, and return a nonzero status on failure.

## `src/geometry.c`

Traced IGRF flux-tube hexahedra in ECEF metres. Apex longitude and apex spherical altitude label a tube; the along-field coordinate is a weighted tracing distance, not geographic latitude.

* [`tangent`](../src/geometry.c#L17): Return the unit magnetic-field tangent in ECEF for field-line integration.
* [`rk4`](../src/geometry.c#L24): Advance one signed physical-distance field-line step with classical RK4 using unit field tangents.
* [`radial`](../src/geometry.c#L32): Return B dot radial-unit-vector; its zero locates a spherical-altitude apex.
* [`tube_apex`](../src/geometry.c#L39): Trace from a given position toward increasing radius and bracket/bisect the radial-field zero; reject escaping or interior trajectories.
* [`half_trace`](../src/geometry.c#L71): Trace one conjugate half tube from its apex to the spherical base boundary, accumulating the altitude-weighted along-field coordinate.
* [`tube_apex_label`](../src/geometry.c#L102): Find the magnetic-apex latitude at specified geographic apex longitude and spherical altitude by a radial-field root solve.
* [`tube_trace`](../src/geometry.c#L120): Trace both conjugate halves and resample ns+1 nodes uniformly in the weighted coordinate; output spans the two base-altitude endpoints.
* [`nid`](../src/geometry.c#L169): Flatten node coordinates with the along-field coordinate contiguous.
* [`tube_cell_id`](../src/geometry.c#L175): Flatten cell coordinates in [ny,nx,ns] storage with k contiguous.
* [`tube_face_id`](../src/geometry.c#L182): Flatten an axis-specific staggered face index, accounting for the extra face in its normal direction.
* [`tube_face_count`](../src/geometry.c#L188): Return the number of faces in one staggered axis family.
* [`jacobian`](../src/geometry.c#L195): Differentiate the trilinear eight-node cell mapping to obtain its three covariant basis vectors.
* [`face_center`](../src/geometry.c#L210): Average the four corners of a staggered face in physical ECEF coordinates.
* [`tube_geometry_allocate`](../src/geometry.c#L226): Allocate nodes/cells/faces for a grid and initialize metric diagnostics; no field tracing is performed here.
* [`tube_geometry`](../src/geometry.c#L241): Trace every node line then construct cell/face metrics over the full grid.
* [`tube_geometry_rows`](../src/geometry.c#L257): Build metrics for a disjoint cell-row range: Jacobian-volume quadrature, transverse covectors, field-aligned tangent, gravity and integrated face flux factors. Reject folded/degenerate geometry.
* [`tube_free`](../src/geometry.c#L363): Free every canonical geometric array and clear the grid object.
* [`chapman`](../src/geometry.c#L374): Evaluate the code's normalized Chapman-like layer exp(1-z-exp(-z)), with z=(height-peak)/scale.
* [`gravity_drift`](../src/geometry.c#L382): Solve the steady perpendicular collisional gravity drift for a species with signed gyrofrequency.
* [`conductivity_cell`](../src/geometry.c#L390): Build O+/fixed-molecular Pedersen/Hall and gravity-current responses, including electrons under quasineutrality.
* [`tube_conductivity`](../src/geometry.c#L419): Populate the explicitly idealized analytic height-dependent collision/conductivity background.
* [`tube_conductivity_table`](../src/geometry.c#L433): Read a strictly increasing height/collision table, require full geometry coverage, log-interpolate collision rates and construct responses.
* [`tube_conductivity_values`](../src/geometry.c#L484): Construct the same conductivity response from externally supplied per-canonical-cell collision frequencies.

## `src/igrf.c`

Independent Schmidt semi-normalized IGRF evaluator. Position is ECEF metres, returned field is ECEF tesla. Geodetic WGS84 height differs from spherical radius minus 6371.2 km.

* [`vadd`](../src/igrf.c#L18): Return component-wise ECEF vector addition.
* [`vsub`](../src/igrf.c#L24): Return component-wise ECEF vector subtraction.
* [`vmul`](../src/igrf.c#L30): Multiply an ECEF vector by a scalar.
* [`vdot`](../src/igrf.c#L36): Return the Euclidean dot product.
* [`vcross`](../src/igrf.c#L42): Return the right-handed Euclidean cross product.
* [`vnorm`](../src/igrf.c#L48): Return the Euclidean vector magnitude.
* [`vunit`](../src/igrf.c#L54): Normalize a nonzero vector; return the zero vector for zero magnitude.
* [`igrf_load`](../src/igrf.c#L62): Read all 195 official IGRF-14 Gauss coefficients, interpolate five-year epochs or apply 2025-2030 secular variation, and retain the requested harmonic degree.
* [`igrf_field`](../src/igrf.c#L118): Evaluate Schmidt semi-normalized harmonics and their derivatives, convert spherical field to ECEF tesla, and use a continuous near-pole evaluation.
* [`spherical_position`](../src/igrf.c#L171): Convert geocentric latitude/longitude in degrees and radius in metres to ECEF metres.
* [`geodetic_position`](../src/igrf.c#L179): Convert WGS84 geodetic latitude/longitude/altitude to ECEF metres.
* [`geodetic_coordinates`](../src/igrf.c#L191): Recover WGS84 geodetic latitude/longitude/altitude from ECEF coordinates by iterative latitude refinement.
* [`field_ned`](../src/igrf.c#L205): Project an ECEF field into local north/east/down components at specified geographic angles.

## `src/integrator.c`

Exponential time differencing for fixed nonnegative collision drag; density uses the matching two-stage explicit update.

* [`exponential_weights`](../src/integrator.c#L12): Evaluate exp(-nu*dt), phi1 and phi2 using stable small-argument series and expm1. Reusing these weights is valid only for the same nu and attempted dt.

## `src/main.c`

Command-line parsing and dispatch. The executable defaults are deliberately distinguished from fully specified example cases; argument validation precedes model construction.

* [`number`](../src/main.c#L22): Parse one finite numeric argument with complete trailing-character/error checks.
* [`integer`](../src/main.c#L33): Parse a nonnegative integer-valued argument in the supported command-line range.
* [`directory`](../src/main.c#L43): Create the output directory on rank zero and reject a directory that already contains run metadata; never overwrite an existing run.
* [`main`](../src/main.c#L68): Initialize MPI/PETSc, verify the thread support needed by OpenMP, parse/validate all CLI controls, load external driver metadata and dispatch to static or AMR evolution.

## `src/map_cli.c`

Geometry inspection program; it traces tubes and reports the Hong Kong target mapping, independently of plasma evolution.

* [`main`](../src/map_cli.c#L18): Parse this executable's interface, validate inputs, run its diagnostic/model entry point, and return a nonzero status on failure.

## `src/neutral_grid.c`

Strict reader for externally supplied frozen neutral/collision grids. Runtime does not evaluate MSIS or an energy equation.

* [`neutral_grid_hash`](../src/neutral_grid.c#L14): Compute the frozen-neutral input file checksum used in physical restart/cache signatures.
* [`neutral_grid_read`](../src/neutral_grid.c#L26): Validate binary magic, dimensions, positions, temperatures, finite neutral/collision data and checksum against canonical geometry; return owned allocated input values.

## `src/night.c`

Local reaction and solar-geometry kernels; no MPI or mutable global state. Reaction density units are m^-3, rates s^-1 or m^3 s^-1, and time seconds.

* [`night_rates`](../src/night.c#L13): Calculate O+-to-NO+/O2+ conversion and dissociative recombination rates from prescribed temperatures and N2/O2; reject unsupported thermal ranges.
* [`convert`](../src/night.c#L38): Apply the cached exact linear conversion operator, splitting lost O+ between molecular products and accumulating converted particles.
* [`recombine`](../src/night.c#L50): Apply the exact scalar solution of dn/dt=-alpha*n*(other+n) for fixed other species, with a stable zero-other limit.
* [`night_react`](../src/night.c#L62): Compose symmetric conversion/NO+/O2+/NO+/conversion operators, retaining nonnegativity and an explicit recombined-electron budget; no photoionization is included.
* [`night_solar_cosine`](../src/night.c#L98): Estimate the solar zenith cosine from elapsed UTC and geographic angles for the declared night-only guard, not precise sunrise or MLT.

## `src/night_probe.c`

Small local chemistry diagnostic used to inspect rate/charge behavior independently of spatial transport.

* [`main`](../src/night_probe.c#L14): Parse this executable's interface, validate inputs, run its diagnostic/model entry point, and return a nonzero status on failure.

## `src/night_transport.c`

Local species mobility, pressure and collisional parallel-flux coefficients. Exact vacuum has a bounded secant continuation rather than a positive density floor.

* [`night_mobility`](../src/night_transport.c#L12): Evaluate bounded Hall/Pedersen mobility factors with hypot-based scaling.
* [`night_log_jump`](../src/night_transport.c#L21): Return log(right)-log(left) for positive states and a bounded secant continuation only at exact vacuum.
* [`B`](../src/night_transport.c#L32): Evaluate the fitted-flux Bernoulli function robustly near zero and at large positive/negative arguments.
* [`night_sg_coeff`](../src/night_transport.c#L47): Compute collisional relaxation and Scharfetter-Gummel-style parallel pressure/gravity coefficients from current electron density.
* [`night_pressure_force`](../src/night_transport.c#L63): Return the parallel pressure acceleration kB/m*(Ti*grad(log ni)+Te*grad(log Ne)).

## `src/output.c`

Uniform-grid metadata, field outputs and checkpoints. Always use saved metadata to decode shapes and units.

* [`open3`](../src/output.c#L14): Open a named uniform-grid output file in the run directory; abort on an I/O failure.
* [`doubles`](../src/output.c#L25): Write exactly n binary double values and reject short writes.
* [`collision_input`](../src/output.c#L43): Copy/hash the effective collision input so the saved run records its actual background policy.
* [`solver3_output`](../src/output.c#L68): Gather uniform conserved states, write fields/metadata/diagnostics and an atomically replaced portable checkpoint.

## `src/pdm.c`

Physical-volume four-cell cubic reconstruction and partial-donor-cell limiting. PDM4 describes the spatial stencil; the coupled time update remains second order.

* [`node_id`](../src/pdm.c#L19): Flatten node coordinates with the along-field coordinate contiguous.
* [`determinant`](../src/pdm.c#L27): Evaluate the signed Jacobian determinant of a trilinear hexahedron at one logical quadrature point.
* [`tube_cell_moments`](../src/pdm.c#L45): Integrate normalized physical-volume moments up to cubic degree using tensor three-point Gaussian quadrature; reject folded cells and inconsistent volumes.
* [`tube_moment_weights`](../src/pdm.c#L87): Solve the four-by-four moment matching system for a face value, using pivoted elimination and polynomial reproduction checks.
* [`tube_pdm_weights`](../src/pdm.c#L135): Prepare four-cell cubic face weights from physical moments for a uniform logical grid; unavailable boundary stencils remain NaN and use a lower-order branch.
* [`tube_pdm_limit`](../src/pdm.c#L189): Clamp the high-order state to donor/downstream bounds, suppress extrema and cap the correction by 2*beta times the upstream difference.

## `src/petsc_potential.c`

Distributed PETSc backend for the uniform-grid reference solver; the AMR solver has its own constrained-mesh assembly.

* [`default_option`](../src/petsc_potential.c#L27): Set a PETSc backend default only when that option has not already been explicitly set.
* [`solver3_petsc_create`](../src/petsc_potential.c#L38): Allocate distributed uniform-grid PETSc matrices/vectors/scatter and configure the reference potential solve.
* [`solver3_petsc_solve`](../src/petsc_potential.c#L78): Assemble current uniform-grid coefficients, apply boundary conditions, solve and verify the true residual before distributing potential.
* [`solver3_petsc_seed`](../src/petsc_potential.c#L152): Warm-start distributed uniform potential DOFs from the current canonical phi.
* [`solver3_petsc_destroy`](../src/petsc_potential.c#L167): Destroy all PETSc objects associated with the uniform-grid backend.

## `src/potential.c`

Small/reference sparse transverse potential solver with bilinear elements and BiCGSTAB/ILU. The production AMR path uses the separate distributed PETSc implementation.

* [`boundary`](../src/potential.c#L14): Classify a potential node or apply its prescribed boundary role; inspect the explicit regional/periodic branches below.
* [`canonical`](../src/potential.c#L21): Map a periodic or physical transverse node to its canonical independent index.
* [`create`](../src/potential.c#L28): Allocate the reference sparse potential graph and solver work arrays for the requested periodic policy.
* [`potential_create`](../src/potential.c#L88): Create a reference transverse operator with nonperiodic longitude.
* [`potential_create_periodic`](../src/potential.c#L94): Create a reference transverse operator with explicit longitude periodicity.
* [`potential_free`](../src/potential.c#L100): Release sparse graph, factorization and iteration arrays.
* [`potential_zero`](../src/potential.c#L114): Clear current matrix/RHS values before assembling a new plasma state.
* [`entry`](../src/potential.c#L121): Locate and accumulate an entry in the preallocated sparse graph.
* [`potential_element`](../src/potential.c#L131): Assemble one bilinear transverse element with both tensor cross terms and gravity-current forcing.
* [`potential_boundary`](../src/potential.c#L159): Apply the reference solver's explicit physical/periodic potential boundary constraints.
* [`apply`](../src/potential.c#L174): Multiply the current sparse potential matrix by an input vector.
* [`dot`](../src/potential.c#L184): Evaluate a vector inner product for the reference iterative solver.
* [`maxabs`](../src/potential.c#L193): Return the infinity norm of a vector.
* [`ilu`](../src/potential.c#L202): Build the reference incomplete-LU preconditioner and reject unusable pivots.
* [`precondition`](../src/potential.c#L226): Apply the stored incomplete-LU triangular factors.
* [`potential_solve`](../src/potential.c#L242): Solve the assembled nonsymmetric system with preconditioned BiCGSTAB and a checked residual criterion.
* [`potential_coefficients`](../src/potential.c#L313): Integrate local conductivity/current responses along field lines for the reference transverse equation.
* [`potential_assemble`](../src/potential.c#L338): Assemble every transverse element and boundary condition from the current integrated coefficients.

## `src/restart.c`

Uniform-grid portable checkpoint reader and physical/numerical input signature. It is distinct from native structure-AMR restart.

* [`read64`](../src/restart.c#L17): Read one little-endian 64-bit unsigned checkpoint word; short reads are fatal.
* [`real64`](../src/restart.c#L29): Decode one checkpoint IEEE double from its little-endian bit pattern.
* [`solver3_file_hash`](../src/restart.c#L38): Calculate an exact input-file FNV identity, with zero for an absent optional input.
* [`array64`](../src/restart.c#L60): Read a portable array of doubles into caller storage.
* [`solver3_checkpoint_parameters`](../src/restart.c#L75): Construct the fixed-order uniform-grid physical/numerical signature for compatibility checks.
* [`solver3_restart`](../src/restart.c#L114): Restore and validate the uniform-grid checkpoint, including original conservation history and imposed-field compatibility.

## `src/transport.c`

Uniform-grid reference solver. Density is normalized by 1e12 m^-3 and momentum is n*u_parallel. This path is retained for static-grid and backend comparisons.

* [`solver3_load_field`](../src/transport.c#L17): Load IGRF, degree-one dipole or the explicit aligned-dipole fixture from the selected magnetic policy.
* [`solver3_initial_density`](../src/transport.c#L30): Evaluate the selected idealized single/double Chapman O+ altitude profile in normalized units.
* [`fail3`](../src/transport.c#L39): Report a fatal rank-local model error and abort MPI when enabled.
* [`alloc3`](../src/transport.c#L49): Allocate zero-initialized reference-solver storage and fail on allocation errors.
* [`wall3`](../src/transport.c#L58): Return a monotonic wall clock for timing; this is independent of model physical time.
* [`solver3_gather`](../src/transport.c#L70): Collect distributed uniform owned state slabs into replicated global arrays when needed.
* [`exchange_state`](../src/transport.c#L84): Exchange uniform-grid density/momentum halos while retaining owned and ghost consistency.
* [`broadcast`](../src/transport.c#L123): Broadcast a potentially large raw array in MPI-count-safe chunks.
* [`solver3_create`](../src/transport.c#L138): Construct uniform geometry, prescribed background, decomposition, reconstruction caches and potential backend.
* [`solver3_free`](../src/transport.c#L315): Release uniform solver fields, geometry, halos, reconstruction weights and potential objects.
* [`solver3_valid`](../src/transport.c#L359): Collectively verify finite, nonnegative uniform density and consistent vacuum momentum.
* [`solver3_mass`](../src/transport.c#L374): Volume-integrate uniform density and sum owned contributions across MPI ranks.
* [`ph`](../src/transport.c#L386): Read a transverse potential node using this solver's boundary/periodic indexing convention.
* [`derivatives`](../src/transport.c#L393): Compute transverse potential differences at a uniform geometric face with its boundary convention.
* [`minmod`](../src/transport.c#L419): Return the sign-consistent slope with smallest magnitude; return zero at an extremum.
* [`bernoulli`](../src/transport.c#L426): Evaluate the Bernoulli fitting function z/(exp(z)-1) with small-argument series and large-argument limits.
* [`value`](../src/transport.c#L440): Read a clamped or periodic scalar state sample.
* [`state`](../src/transport.c#L463): Reconstruct a uniform face state using the selected slope policy.
* [`exb_state`](../src/transport.c#L475): Reconstruct a bounded transverse density face value with physical-volume PDM4 or the configured lower-order stencil.
* [`velocity_value`](../src/transport.c#L498): Return parallel velocity with an exact-vacuum zero branch.
* [`velocity_state`](../src/transport.c#L507): Construct the bounded parallel-velocity face state used for momentum transport.
* [`line_length`](../src/transport.c#L519): Read the physical length of an along-field segment.
* [`log_difference`](../src/transport.c#L527): Compute a positive-state logarithmic pressure jump or its exact-vacuum secant extension.
* [`pressure_log_gradient`](../src/transport.c#L534): Construct a nonuniform parallel logarithmic density gradient for the uniform isothermal pressure force.
* [`solver3_rhs`](../src/transport.c#L555): Exchange current state, integrate conductivity, solve the potential, and evaluate the uniform transport/source RHS.
* [`solver3_transport`](../src/transport.c#L620): Evaluate conservative uniform particle/momentum face fluxes, gravity/pressure sources, boundary budgets and CFL rates.

## Inline ownership helpers

`night_start(a)` is the global first-owned canonical offset. `night_owned(a)` is its owned-cell count without ghosts. Both are directly commented in `src/amr3_night_internal.h`.
