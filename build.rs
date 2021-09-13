use std::process::Command;
use std::env;
use std::path::Path;

fn main() {
    let target_arch = env::var("CARGO_CFG_TARGET_ARCH").unwrap();
    let out_dir = env::var("OUT_DIR").unwrap();

    // Note that there are a number of downsides to this approach, the comments
    // below detail how to improve the portability of these commands.
    Command::new("riscv64-unknown-elf-gcc").args(&["-c", "src/asm/gfp_rvv.s", "-o", "gfp_rvv.o"])
        .arg(&format!("{}/gfp_rvv.o", out_dir))
        .status().unwrap();
    println!("cargo:rustc-link-search=native={}", out_dir);
    println!("cargo:rerun-if-changed=src/asm/gfp_rvv.s");

    // Command::new("ar").args(&["crus", "libhello.a", "hello.o"])
    //     .current_dir(&Path::new(&out_dir))
    //     .status().unwrap();
    // println!("cargo:rustc-link-lib=static=gfp_rvv");

    // cc::Build::new()
    //     .file("src/asm/gfp_rvv.s")
    //     .compile("gfp_rvv");
}
