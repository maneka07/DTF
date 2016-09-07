#include "pfarb_buf_io.h"
#include "pfarb_file_buffer.h"
#include <assert.h>
#include "pfarb.h"
#include "pfarb_mem.h"

static MPI_Offset to_1d_index(int ndim, MPI_Offset *shape, MPI_Offset *coord)
{
    int dimid = 0;
    MPI_Offset idx = coord[dimid];

    while(dimid != ndim - 1){
       idx = idx*shape[dimid+1] + coord[dimid+1];
       dimid++;
    }
    return idx;
}

MPI_Offset read_write_var(const char *filename, int varid, const MPI_Offset *start,
                          const MPI_Offset *count, const MPI_Offset *stride, const MPI_Offset *imap,
                          MPI_Datatype dtype, void *buf, int rw_flag)
{
    MPI_Offset el_offt, byte_offt, lead_el_idx, el_cnt;
    int dimid;
    int el_sz;
    MPI_Offset data_sz;
    MPI_Offset buf_offt = 0, tmp, ret;

    file_buffer_t *fbuf = find_file_buffer(gl_filebuf_list, filename);
    assert(buf!=NULL);

    if(rw_flag == FARB_READ)
         assert(fbuf->is_ready);

    if(imap != NULL){
        FARB_DBG(VERBOSE_ERROR_LEVEL, "FARB Error: writing mapped vars is not impelemented yet. Ignore.");
        return 0;
    }

    if(stride != NULL){
        FARB_DBG(VERBOSE_ERROR_LEVEL, "FARB Error: writing vars at a stride is not impelemented yet. Ignore.");
        return 0;
    }

    farb_var_t *var = find_var(fbuf->vars, varid);
    if(var == NULL){
        FARB_DBG(VERBOSE_ERROR_LEVEL, "FARB Error: could not find var with id %d", varid);
        return 0;
    }

    MPI_Type_size(dtype, &el_sz);
    assert(el_sz > 0);
    if(var->ndims <= 1){ /*scalar or 1d array*/
        el_offt = *start;
        byte_offt = el_offt*(MPI_Offset)el_sz;
        /*Data size to write*/
        data_sz = (*count)*el_sz;

        if(rw_flag == FARB_READ){
            FARB_DBG(VERBOSE_DBG_LEVEL, "Var %d: will read %llu elems of sz %d from element offt %llu ", var->id, *count, el_sz, el_offt);
            ret = mem_read(var, byte_offt, data_sz, buf);
        } else {
            FARB_DBG(VERBOSE_DBG_LEVEL, "Var %d: will write %llu elems of sz %d to element offt %llu ", var->id, *count, el_sz, el_offt);
            ret = mem_write(var, byte_offt, data_sz, buf);
        }

        if(ret != data_sz)
            FARB_DBG(VERBOSE_DBG_LEVEL, "FARB Warning: Meant to read/write %llu bytes but actual val is %llu", data_sz, ret);

    } else { /*multi-dimentional array*/
        MPI_Offset *start1 = malloc(var->ndims*sizeof(MPI_Offset));
        assert(start1 != NULL);
        memcpy(start1, start, var->ndims*sizeof(MPI_Offset));
        FARB_DBG(VERBOSE_DBG_LEVEL, "start el %lld, until %lld", start1[0], start1[0] + count[0]);
        for(lead_el_idx = start[0]; lead_el_idx < start[0] + count[0]; lead_el_idx++){
            /*Calculate element offset within buffer nodes*/
            /*for an array stored in row-major form, elements with a fixed row coordinate
            (e.g. coord X in XYZ dimensions) will be written contiguously in a 1d array.
            So we just need to find an element offset for a given X and then write Y*Z
            contiguous elements. */
            start1[0] = lead_el_idx;
            el_offt = to_1d_index(var->ndims, var->shape, start1);
            byte_offt = el_offt*el_sz;

            /*how many elements to write contiguously?*/
            dimid = 1;
            el_cnt = count[dimid];
            while(dimid != var->ndims - 1){
                el_cnt = el_cnt*count[dimid+1];
                dimid++;
            }

            data_sz = el_cnt*el_sz;

            if(rw_flag == FARB_READ){
                FARB_DBG(VERBOSE_DBG_LEVEL, "Var %d: will read %llu elems of sz %d from element offt %llu ", var->id, el_cnt, el_sz, el_offt);
                tmp = mem_read(var, byte_offt, data_sz, (void*)(buf+buf_offt));
            } else { /*FARB_WRITE*/
                FARB_DBG(VERBOSE_DBG_LEVEL, "Var %d: will write %llu elems of sz %d to element offt %llu ", var->id, el_cnt, el_sz, el_offt);
                tmp = mem_write(var, byte_offt, data_sz, (void*)(buf+buf_offt));
            }

            if(tmp != data_sz)
                FARB_DBG(VERBOSE_DBG_LEVEL, "FARB Warning: Meant to read/write %llu bytes but actual val is %llu", data_sz, tmp);
            buf_offt+=tmp;
        }
        free(start1);
        ret = tmp;
    }

    return ret;
}
