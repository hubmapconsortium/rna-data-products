FROM hubmap/azimuth-annotate:1.8
COPY bin/azimuth.R /opt
COPY bin/mtx_to_seurat.R /opt
RUN chmod -R +x /opt/azimuth.R
RUN chmod -R +x /opt/mtx_to_seurat.R 
CMD ["/bin/bash"]