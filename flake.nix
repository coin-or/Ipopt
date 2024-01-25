{
  description = "IPOPT";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
    mumps = {
      url = "github:dzmitry-lahoda-forks/mumps/2adfc01cd0306b30379886cfc0746e5135488a92";
    };
  };

  outputs = inputs@{ flake-parts, ... }:
    flake-parts.lib.mkFlake { inherit inputs; } {
      systems = [ "x86_64-linux" "aarch64-darwin" "x86_64-darwin" ];
      perSystem = { config, self', inputs', pkgs, system, ... }:
        {
          packages = with pkgs;
            rec {
              ipopt = pkgs.stdenv.mkDerivation {
                name = "ipopt";
                version = "3.14.13";
                src = ./.;
                CXXDEFS = [
                  "-D HAVE_RAND"
                  "-D HAVE_CSTRING"
                  "-D HAVE_CSTDIO"
                ];

                configureFlags = [
                  "--with-mumps-cflags=-I${inputs'.mumps.packages.mumps}/include"
                  "--with-mumps-lflags=libdmumps"

                ];

                nativeBuildInputs = [ pkg-config gfortran ];
                buildInputs = [ blas lapack openmpi inputs'.mumps.packages.mumps ];

                enableParallelBuilding = true;
              };
              default = ipopt;

            };
        };
    };
}
