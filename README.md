# cfhtProc2caom2
An application to generate CAOM2 Observations from processed CFHT (NGVS and MegaPipe) FITS files.

# How To Run cfhtProc2caom2

In an empty directory (the 'working directory'), on a machine with Docker installed:

1. In the master branch of this repository, find the scripts directory, and copy the file `cfht_proc_run.sh` to the working directory. e.g.:

  ```
  wget https://raw.github.com/opencadc-metadata-curation/cfhtProc2caom2/master/scripts/cfht_proc_run.sh
  ```

2. Ensure the script is executable:

```
chmod +x cfht_proc_run.sh
```

3. To run the application:

```
./cfht_proc_run.sh
```
Running this will ingest the records in the file `todo.txt`.
