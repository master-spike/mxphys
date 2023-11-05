# MxPhys
MxPhys is a 2d rigid body physics engine.


## Building

To build the demo, clone the git repository into somewhere on your local machine:

```git clone git@github.com:master-spike/mxpys```

```cd mxphys```

Then follow the platform specific instructions.

### Linux

Make a directory called `build` inside where you cloned MxPhys to.

```mkdir build```

Then do `cd build` and `cmake ../` and `cmake --build .` to build mxphys. If this doesn't work, you may need to do

`apt-get install libsdl2-dev`

to get the sdl2 libraries on your system.

### Windows (Visual Studio)

Ensure you have git and cmake installed on your system, and that you use Powershell as your terminal.

Download SDL2 from https://github.com/libsdl-org/SDL/releases/tag/release-2.28.5, choosing the version which ends in `VC`, and extract it to some location.
Note the filepath to the `cmake/` directory.

Do the same steps as for linux, except when doing `cmake ../`, you must pass `-DSDL2_DIR=path/to/sdl/cmake` so that cmake can find SDL.

After doing `cmake --build .`, in order to run the `mxphysdemo.exe` executable, you must copy `SDL2.dll` to the same location.

## Include MxPhys in your own project

If you project uses git, you can add this project as a git submodule to `thirdparty/mxphys` for example, and in your top level `CMakeLists.txt`,
add the lines `add_subdirectory(path/to/mxphys)` and for each target which includes mxphys, `target_include_directories(path/to/mxphys/include)`
and add `MxPhysLib` to your `target_link_libraries`.

## Contributing

Contributions in the form of pull requests are generally welcome. If you want a feature to be added, or encounter an issue, submit an issue
stating your request/issue. To contribute code, fork the repository and create feature branches in your fork for any features you wish to
contribute. Then Pull Request to the main branch here.


