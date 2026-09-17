# Container image for RSeQC-rust (PORTING_PLAN Step 10: "container
# images as appropriate"). Multi-stage build: compile in a full Rust
# image, ship only the resulting binaries in a minimal runtime image.
#
# NOTE: this Dockerfile has not been built/tested in this development
# sandbox (no `docker` available here) -- verified instead by checking
# each base image tag actually exists (Docker Hub API) and that the
# workspace has no C-toolchain build dependencies (see
# .github/workflows/ci.yml's own note on the same point). Treat as a
# reasonable first draft, not a confirmed-working artifact, until it's
# actually built once.
#
# Build:  docker build -t rseqc-rust .
# Run:    docker run --rm -v "$PWD:/data" rseqc-rust bam_stat.py -i /data/sample.bam

FROM rust:1.98.1-slim-bookworm AS builder
WORKDIR /build
COPY . .
# Stage just the 33 compiled binaries (not target/release/'s deps/,
# build/, *.d files, etc.) plus their original-upstream-name aliases
# into a clean directory, same "find the top-level executable files"
# pattern scripts/build-release-archive.sh already uses -- keeps the
# runtime image to just what it needs.
RUN cargo build --workspace --release --locked \
    && mkdir -p /build/dist-bin \
    && find target/release -maxdepth 1 -type f -executable -exec cp {} /build/dist-bin/ \; \
    && ./scripts/install-aliases.sh /build/dist-bin --copy

FROM debian:bookworm-slim
# rseqc-formats' noodles-* dependencies do their own decompression in
# pure Rust (no libz/libbz2 system package needed -- see CI's own note
# on this); nothing beyond libc is required here.
COPY --from=builder /build/dist-bin/ /usr/local/bin/
ENTRYPOINT []
CMD ["bam_stat.py", "--help"]
