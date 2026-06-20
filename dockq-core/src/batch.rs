//! Batch driver (task #6, integration owner). Two first-class shapes, both parallel in
//! Rust (Rayon): one-native-vs-many-models, and arbitrary (model, native) pair lists,
//! plus a directory-scan convenience. Parsing + scoring of each job runs entirely in Rust.
//! Errors per job are reported explicitly (no silent skips).
