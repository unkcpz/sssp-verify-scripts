## Import a new library as a group

I assume the new library named `x_lib`.

- create a folder `x_lib` and put the UPF file inside it.
- update the `src/sssp_verify_scripts/scripts/import_upf_group.py` with a new `x_lib` library item.
- run `import_upf_group --base-folder ~/project/sssp-project/sssp-verify-scripts --lib-name x_lib`
