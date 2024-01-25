{
  description = "IPOPT";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
    mumps = {
      url = "github:dzmitry-lahoda-forks/mumps/e6df3ef20e00e6be5c5e0a859ac517bbfd28827b";
    };
  };

  outputs = inputs@{ flake-parts, ... }:
    flake-parts.lib.mkFlake { inherit inputs; } {
      systems = [ "x86_64-linux" "aarch64-darwin" "x86_64-darwin" ];
      perSystem = { config, self', inputs', pkgs, system, ... }:
        {
          packages = with pkgs;
            rec {
              ipopt-mumps-seq = pkgs.stdenv.mkDerivation {
                name = "ipopt";
                version = "3.14.13";
                src = ./.;
                CXXDEFS = [
                  "-D HAVE_RAND"
                  "-D HAVE_CSTRING"
                  "-D HAVE_CSTDIO"
                ];

                configureFlags = [
                  "--with-mumps-cflags=-I${inputs'.mumps.packages.mumps-32-seq}/include"
                  "--with-mumps-lflags=libdmumps"
                  "--disable-mpiinit"

                ];

                nativeBuildInputs = [ pkg-config gfortran ];
                buildInputs = [ blas lapack inputs'.mumps.packages.mumps-32-seq ];

                enableParallelBuilding = true;
              };
              default = ipopt-mumps-seq;

            };
        };
    };
}
