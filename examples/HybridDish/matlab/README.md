MATLAB readout for HybridDish
=============================

BSim writes semicolon-delimited CSVs, the same convention as reservoir_new.
Each run folder also has matlab_meta.txt with reshape instructions for one
sample row (20x10 state maps, 4x2 count maps).

Typical session, from the HybridDish folder:

    addpath('matlab')
    S = ridge_narma10('output/driven_seed101')

    % Controls / diagnostics
    Sb = ridge_narma10('output/brownian_seed101', 'input/narma10_target.csv', ...
                       'Family', 'brownian');
    Sf = ridge_narma10('output/driven_seed101', 'input/narma10_target.csv', ...
                       'Family', 'field');

    W = classify_waveform('output/waveform_driven_seed101');

Wringe pack 2 (waveform F1 / confusion on frozen waveform voxels;
does not run BSim). Lambda is highest macro OVR val AUC, not NRMSE:

    R = report_waveform_pack
    R = report_waveform_pack('Only', "driven_seed101")

See examples/BSimReservoirPlanWaveform/results/WRINGE_WAVEFORM_PACK.md.

Wringe pack 1 (NARMA-10 reporting layer on frozen Narma10b voxels;
does not run BSim):

    R = report_narma_pack
    R = report_narma_pack('Only', "driven_seed222")

    B  = baselines(target.y_next)           % intercept + persistence
    MC = memory_capacity(X, target.u)       % linear MC, k=1..20
    H  = horizon_narma(X, target.y_next, 5) % teacher-forced h-step

See examples/BSimReservoirPlanNarma10b/results/WRINGE_NARMA_PACK.md.

Wringe pack 3 (Mackey–Glass one-pager on frozen BenchA voxels;
does not run BSim). Teacher-forced one-step x[n+1], not autonomous MG:

    R = report_mg_pack
    R = report_mg_pack('Only', "driven_seed101")

    B  = baselines(target.x_next, 'Label', 'Mackey-Glass')
    Sf = ridge_closed(Xf, target.x_next)   % field-only diagnostic

See examples/BSimReservoirPlanBenchA/results/WRINGE_MG_PACK.md.

Wringe pack 4 (kernel-rank diagnostic on frozen Narma10b voxels;
does not run BSim). 20×10 maps with S = N = 200. Official 408-D is a
ceiling, not a KR claim. Not CHARC. Not a gate:

    R = report_kr_pack
    R = report_kr_pack('Only', "driven_seed222")

    [Xl, ~] = window_features(T, 'lum');
    Kl = kernel_rank(Xl, 'Name', 'Lum_Mean_*');

See examples/BSimReservoirPlanNarma10b/results/WRINGE_KR_PACK.md.

Official biology features (408):
  Receiver_R_*     window mean, 20 x 10
  Lum_Mean_*       window mean, 20 x 10
  Input_Driven_Death_*   last sample of the window, 4 x 2

Voxel AHL_uM_* is a field-only diagnostic. It is not the official readout.

Ridge split (frozen):
  washout 0..39, train 40..149, test 150..199
  lambda chosen on train windows 128..149 only
