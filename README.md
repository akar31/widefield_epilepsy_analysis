# Widefield_epilepsy_analysis
### Azra Karatan
### 06/23/2026

## Directory Layout:
```
epilepsy_01/
|
|  |
|  |__/Analyze_one_WF_run.m
|  |__/applyCortexMask.m
|  |__/brain_mask.mat
|  |__/get_signal.m
|  |__/seizure_analysis1.m
|  |

```

NOTES:
 - Use seizure_analysis1.m to create/apply cortex mask. Mask already exists, no need to overwrite unless changes needed.
     * complete code
       
 - Use get_signal.m to compute the DeltaF/F signal from the pixel fluorescence values.
     * currently runs for 8 frames for testing
     * f0 calculation will be fixed to compute movie_blue, movie_green separately. 


## OS and C compiler
OS: MacOS Sequoia 15.6.1
C compiler: Apple clang version 17.0.0 (clang-1700.6.3.2)


