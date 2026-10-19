### CPU time per Newton step

Collected 2026-10-19: 16 of 18 cases completed. Block-ILU(0) at k=7 and k=8
failed during DG with OUT_OF_MEMORY at a 120 GiB allocation; HDG was not run.
Failed-process RSS below is not a successful-run memory requirement.

CPU: Intel(R) Xeon(R) Gold 6242 CPU @ 2.80GHz (node1601), one thread.
Despite the historical table title, mean_newton_s measures elapsed wall time
per accepted Newton update, not process CPU time. Process CPU and wall times
are recorded separately in MEMORY_ESTIMATE.md and each case's resources.csv.

Configuration: 256x256, CFL=5.0, SSP-DIRK3, four stages, three fixed
Newton updates per stage (12 measured updates), matrix-free partial
assembly, restarted weighted right GMRES (restart dimension 20, maximum 1500 Krylov iterations), and one process/thread.

Peak RSS is measured per solver process over its complete lifetime,
including setup; missing measurements are shown as —, not estimates.

| k | preconditioner | CFL | DG CPU time per Newton step (s) | HDG CPU time per Newton step (s) | DG/HDG | DG peak RSS (GiB) | HDG peak RSS (GiB) | node | job_id |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0 | block-ilu0 | 5.0 | 0.63753157566396523 | 2.0773129234904943 | 0.30690204083105854 | 0.2716484069824219 | 0.569091796875 | node1601.oscar.ccv.brown.edu | 7583445 |
| 0 | block-jacobi | 5.0 | 1.2944837657507728 | 3.4374422644117448 | 0.3765834205137697 | 0.2038726806640625 | 0.3782844543457031 | node1601.oscar.ccv.brown.edu | 7583444 |
| 1 | block-ilu0 | 5.0 | 3.5944653986586501 | 5.5809946404750947 | 0.6440546229144315 | 1.2134208679199219 | 1.4156227111816406 | node1601.oscar.ccv.brown.edu | 7583443 |
| 1 | block-jacobi | 5.0 | 6.0816388836003581 | 7.5270784365879848 | 0.8079680496005509 | 0.6711692810058594 | 0.9375953674316406 | node1601.oscar.ccv.brown.edu | 7583442 |
| 2 | block-ilu0 | 5.0 | 15.661668416753772 | 14.593375791411079 | 1.0732039413369616 | 4.321250915527344 | 2.917133331298828 | node1601.oscar.ccv.brown.edu | 7583441 |
| 2 | block-jacobi | 5.0 | 27.856269659506665 | 24.902720978916022 | 1.118603452333232 | 1.7382240295410156 | 1.9668388366699219 | node1601.oscar.ccv.brown.edu | 7583440 |
| 3 | block-ilu0 | 5.0 | 40.146268260489116 | 24.858708220243546 | 1.614978055367987 | 11.949943542480469 | 5.264698028564453 | node1601.oscar.ccv.brown.edu | 7583439 |
| 3 | block-jacobi | 5.0 | 59.157091871153646 | 38.210232769404904 | 1.5482002485606678 | 3.8802223205566406 | 3.654285430908203 | node1601.oscar.ccv.brown.edu | 7583438 |
| 4 | block-ilu0 | 5.0 | 109.98157531642967 | 45.723570625073414 | 2.4053584138969044 | 27.354736328125 | 8.634872436523438 | node1601.oscar.ccv.brown.edu | 7583437 |
| 4 | block-jacobi | 5.0 | 144.84621904102275 | 82.978139214334121 | 1.7455949291280466 | 7.740577697753906 | 6.175926208496094 | node1601.oscar.ccv.brown.edu | 7583436 |
| 5 | block-ilu0 | 5.0 | 218.36511683250339 | 72.284765028669895 | 3.0209009705695866 | 54.7783203125 | 13.223461151123047 | node1601.oscar.ccv.brown.edu | 7583435 |
| 5 | block-jacobi | 5.0 | 251.45108258575462 | 125.85146847617034 | 1.9979987967591037 | 14.17373275756836 | 9.731624603271484 | node1601.oscar.ccv.brown.edu | 7583434 |
| 6 | block-ilu0 | 5.0 | 445.10608194175228 | 116.70809940740641 | 3.8138405492147482 | 99.36463165283203 | 19.210681915283203 | node1601.oscar.ccv.brown.edu | 7583433 |
| 6 | block-jacobi | 5.0 | 452.2425671913295 | 212.96592885597963 | 2.123544219588116 | 24.203170776367188 | 14.49472427368164 | node1601.oscar.ccv.brown.edu | 7583432 |
| 7 | block-jacobi | 5.0 | 794.52435649201561 | 353.14934146215506 | 2.2498253945552382 | 39.05031967163086 | 20.670978546142578 | node1601.oscar.ccv.brown.edu | 7583430 |
| 8 | block-jacobi | 5.0 | 1375.5149045302605 | 474.81393736684703 | 2.8969556204655404 | 60.13582992553711 | 28.433704376220703 | node1601.oscar.ccv.brown.edu | 7583428 |

### Measured process peak RSS

Failed processes are retained here; their RSS is not a completed-run estimate.

| k | preconditioner | form | status | peak RSS (GiB) | wall (s) | node | job ID |
|---:|---|---|---|---:|---:|---|---:|
| 0 | block-ilu0 | dg | completed | 0.2716484069824219 | 10.94 | node1601.oscar.ccv.brown.edu | 7583445 |
| 0 | block-ilu0 | hdg | completed | 0.569091796875 | 28.53 | node1601.oscar.ccv.brown.edu | 7583445 |
| 0 | block-jacobi | dg | completed | 0.2038726806640625 | 18.73 | node1601.oscar.ccv.brown.edu | 7583444 |
| 0 | block-jacobi | hdg | completed | 0.3782844543457031 | 44.61 | node1601.oscar.ccv.brown.edu | 7583444 |
| 1 | block-ilu0 | dg | completed | 1.2134208679199219 | 47.55 | node1601.oscar.ccv.brown.edu | 7583443 |
| 1 | block-ilu0 | hdg | completed | 1.4156227111816406 | 71.71 | node1601.oscar.ccv.brown.edu | 7583443 |
| 1 | block-jacobi | dg | completed | 0.6711692810058594 | 77.05 | node1601.oscar.ccv.brown.edu | 7583442 |
| 1 | block-jacobi | hdg | completed | 0.9375953674316406 | 94.65 | node1601.oscar.ccv.brown.edu | 7583442 |
| 2 | block-ilu0 | dg | completed | 4.321250915527344 | 194.77 | node1601.oscar.ccv.brown.edu | 7583441 |
| 2 | block-ilu0 | hdg | completed | 2.917133331298828 | 181.33 | node1601.oscar.ccv.brown.edu | 7583441 |
| 2 | block-jacobi | dg | completed | 1.7382240295410156 | 339.41 | node1601.oscar.ccv.brown.edu | 7583440 |
| 2 | block-jacobi | hdg | completed | 1.9668388366699219 | 304.31 | node1601.oscar.ccv.brown.edu | 7583440 |
| 3 | block-ilu0 | dg | completed | 11.949943542480469 | 494.4 | node1601.oscar.ccv.brown.edu | 7583439 |
| 3 | block-ilu0 | hdg | completed | 5.264698028564453 | 306.78 | node1601.oscar.ccv.brown.edu | 7583439 |
| 3 | block-jacobi | dg | completed | 3.8802223205566406 | 716.92 | node1601.oscar.ccv.brown.edu | 7583438 |
| 3 | block-jacobi | hdg | completed | 3.654285430908203 | 465.88 | node1601.oscar.ccv.brown.edu | 7583438 |
| 4 | block-ilu0 | dg | completed | 27.354736328125 | 1342.49 | node1601.oscar.ccv.brown.edu | 7583437 |
| 4 | block-ilu0 | hdg | completed | 8.634872436523438 | 560.01 | node1601.oscar.ccv.brown.edu | 7583437 |
| 4 | block-jacobi | dg | completed | 7.740577697753906 | 1748.07 | node1601.oscar.ccv.brown.edu | 7583436 |
| 4 | block-jacobi | hdg | completed | 6.175926208496094 | 1005.32 | node1601.oscar.ccv.brown.edu | 7583436 |
| 5 | block-ilu0 | dg | completed | 54.7783203125 | 2660.91 | node1601.oscar.ccv.brown.edu | 7583435 |
| 5 | block-ilu0 | hdg | completed | 13.223461151123047 | 882.3 | node1601.oscar.ccv.brown.edu | 7583435 |
| 5 | block-jacobi | dg | completed | 14.17373275756836 | 3031.77 | node1601.oscar.ccv.brown.edu | 7583434 |
| 5 | block-jacobi | hdg | completed | 9.731624603271484 | 1522.74 | node1601.oscar.ccv.brown.edu | 7583434 |
| 6 | block-ilu0 | dg | completed | 99.36463165283203 | 5409.98 | node1601.oscar.ccv.brown.edu | 7583433 |
| 6 | block-ilu0 | hdg | completed | 19.210681915283203 | 1419.94 | node1601.oscar.ccv.brown.edu | 7583433 |
| 6 | block-jacobi | dg | completed | 24.203170776367188 | 5448.5 | node1601.oscar.ccv.brown.edu | 7583432 |
| 6 | block-jacobi | hdg | completed | 14.49472427368164 | 2572.16 | node1601.oscar.ccv.brown.edu | 7583432 |
| 7 | block-ilu0 | dg | failed | 119.74918746948242 | 87.51 | node1601.oscar.ccv.brown.edu | 7583431 |
| 7 | block-jacobi | dg | completed | 39.05031967163086 | 9566.69 | node1601.oscar.ccv.brown.edu | 7583430 |
| 7 | block-jacobi | hdg | completed | 20.670978546142578 | 4259.16 | node1601.oscar.ccv.brown.edu | 7583430 |
| 8 | block-ilu0 | dg | failed | 119.7475700378418 | 89.5 | node1601.oscar.ccv.brown.edu | 7583429 |
| 8 | block-jacobi | dg | completed | 60.13582992553711 | 16552.07 | node1601.oscar.ccv.brown.edu | 7583428 |
| 8 | block-jacobi | hdg | completed | 28.433704376220703 | 5723.66 | node1601.oscar.ccv.brown.edu | 7583428 |
