# DLL Interface Reference

Analytic-ICP currently uses two dynamically loaded DLL layers:

```text
Test.exe
  └─ SmoothAdjustment.dll
       └─ Analytic_ICP.dll
```

The exported functions use `extern "C" __declspec(dllexport)` and are resolved with Windows `GetProcAddress`. The checked-in implementation uses the default Windows C/C++ calling convention and does not yet provide a versioned ABI.

## Core Interface: `Analytic_ICP.dll`

The core symbols are defined in `Analytic_ICP/dllmain.cpp`.

### 2D registration

```cpp
double feature_transf2d(
    double (*point_set4regist)[2],
    int num4regist,
    double (*aim_point_set)[2],
    int aim_num,
    int iterCount);
```

Transforms `point_set4regist` in place toward `aim_point_set` using 2D Analytic-ICP.

### 3D registration

```cpp
double feature_transf3d(
    double (*point_set4regist)[3],
    int num4regist,
    double (*aim_point_set)[3],
    int aim_num,
    int iterCount);
```

Transforms `point_set4regist` in place toward `aim_point_set` using 3D Analytic-ICP.

### Synthetic 2D analytic deformation

```cpp
void randomly_2d_map(
    double (*target)[2],
    double (*moved)[2],
    int num,
    int deg);
```

Writes a randomly generated 2D analytic deformation of `target` into `moved`.

### Synthetic 3D analytic deformation

```cpp
void randomly_3d_map(
    double (*target)[3],
    double (*moved)[3],
    int num,
    int deg);
```

Writes a randomly generated 3D analytic deformation of `target` into `moved`.

## File Wrapper: `SmoothAdjustment.dll`

The wrapper symbols are defined in `SmoothAdjustment/dllmain.cpp`.

### Initialization

```cpp
bool aicp_initialize(TCHAR* dllPath);
```

Loads `Analytic_ICP.dll` and resolves:

```text
feature_transf2d
feature_transf3d
randomly_2d_map
randomly_3d_map
```

It must succeed before calling the remaining wrapper functions.

### Release

```cpp
void aicp_release();
```

Releases the loaded core DLL and clears the wrapper's global function pointers.

### File-based 2D registration

```cpp
double p2d_regist(
    const char* fixed_path,
    const char* moving_path,
    const char* moved_output);
```

Reads 2D fixed and moving point sets, registers the moving set, writes the transformed moving points to `moved_output`, and returns the scalar value returned by the core algorithm.

### File-based 3D registration

```cpp
double p3d_regist(
    const char* fixed_path,
    const char* moving_path,
    const char* moved_output);
```

Reads 3D fixed and moving point sets, registers the moving set, writes the transformed moving points to `moved_output`, and returns the scalar value returned by the core algorithm.

### File-based synthetic deformation

```cpp
void analytic_2d_map(
    const char* input_path,
    const char* deformed_output,
    int deg);

void analytic_3d_map(
    const char* input_path,
    const char* deformed_output,
    int deg);
```

These wrappers read a point set, invoke the corresponding random analytic deformation function, and write the deformed points.

## Required Call Sequence

A direct wrapper client should:

1. load `SmoothAdjustment.dll`;
2. resolve the required wrapper symbols;
3. call `aicp_initialize` with the path to `Analytic_ICP.dll`;
4. call one or more registration or deformation functions;
5. call `aicp_release`;
6. release `SmoothAdjustment.dll`.

`Test/Test.cpp` is the reference client for this sequence.

## File and Memory Assumptions

- Input files have no header and contain one 2D or 3D point per row.
- The current parser accepts comma or whitespace separators.
- The wrapper allocates fixed arrays for up to 409,600 points.
- The parser does not perform an explicit bounds check before writing into those arrays; inputs must remain below that implementation limit.
- The core array interfaces are mutable. Callers should not assume that arrays passed through the low-level API remain unchanged.
- The libicp-derived matching code requires at least five model points.

## Error Handling

The wrapper returns `-1.0` for several invalid-pointer, initialization, or empty-input conditions. Other I/O failures are not represented by a uniform status type:

- input opening can throw `std::runtime_error`;
- output opening can terminate the process with `exit(EXIT_FAILURE)`;
- the console frontend reports Windows DLL loading errors through `GetLastError`.

Production callers should therefore isolate the prototype boundary and validate all files before invocation.

## Thread Safety and ABI Status

`SmoothAdjustment.dll` stores the core DLL handle and resolved function pointers in global variables. The initialization/release layer is not designed for concurrent independent sessions.

The API exposes compiler-specific Windows types and raw C-style arrays. There is no semantic-versioned binary compatibility guarantee in the current research snapshot.
