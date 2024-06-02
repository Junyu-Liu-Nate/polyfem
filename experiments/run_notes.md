# Run Notes

## Command line

- **coarse-golf-ball-3D-P3**
    ```
    ./PolyFEM_bin -j ../experiments/high-order-ipc-data/scripts/higher-order/coarse-golf-ball-3D-P3.json -o ../experiments/output/contact-3D-higher-order/coarse-golf-ball-3D-P3/
    ```
- **hex-spline**
    ```
    ./PolyFEM_bin -j ../experiments/json/iga-ipc/2d-quad-spline.json -o ../experiments/output/contact-2D-IGA/2d-quad-spline/
    ```

## Results

- **2-cubes-falling-spline**
    ```
    2024-05-29 12:08:27.417] [polyfem] [info] Building collision mesh...
    [2024-05-29 12:08:27.417] [polyfem] [debug] Building collision proxy with max edge length=1.44 ...
    [2024-05-29 12:08:27.417] [polyfem] [error] build_collision_proxy() is only implemented for tetrahedra!
    libc++abi: terminating due to uncaught exception of type std::runtime_error: build_collision_proxy() is only implemented for tetrahedra!
    Abort trap: 6
    ```

- **coarse-golf-ball-3D-P3**
    Success

- **hex-spline**
    Segmentation fault