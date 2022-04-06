#pragma once
#include <iostream>
#include <algorithm>
#include <string>

#include <stdlib.h>    // free
#include "zstd.h"

/* #include <sstream>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/zstd.hpp> */

namespace utility{
    std::pair<void*,size_t> compress(void* input_buffer, size_t input_buffer_size){
        std::cout << "COMPRESSING" << std::endl;
        size_t output_buffer_size_bound = ZSTD_compressBound(input_buffer_size);
        void* const output_buffer = malloc(output_buffer_size_bound);
        std::cout << "output_buffer_size_bound: " << output_buffer_size_bound << std::endl;
        std::cout << "output_buffer before compression: " << (uint64_t) output_buffer << std::endl;
        size_t output_size = ZSTD_compress(output_buffer, output_buffer_size_bound, input_buffer, input_buffer_size, 1);
        std::cout << "output_buffer after compression: " << (uint64_t) output_buffer << std::endl;
        std::cout << "output_size: " << output_size << std::endl;
        if (ZSTD_isError(output_size)){
            std::cerr << "Compression Error: - " << ZSTD_getErrorName(output_size) << std::endl;
            exit(-1);
        }
        return std::pair<void*,size_t>(output_buffer,output_size);
    }

    size_t decompress(void* input_buffer, void* output_buffer, size_t compressed_size){
        std::cout << "DECOMPRESSING" << std::endl;
        
        uint64_t rSize = ZSTD_getFrameContentSize(input_buffer, compressed_size);
        if (rSize == ZSTD_CONTENTSIZE_ERROR){
            std::cerr << "Decompression Error: File was not compressed by ZSTD" << std::endl;
            return -1;
        }

        if(rSize == ZSTD_CONTENTSIZE_UNKNOWN){
            std::cerr << "Decompression Error: File unable to get content size" << std::endl;
            return -1;
        }
        

        return ZSTD_decompress(output_buffer, rSize, input_buffer, compressed_size);
        /* int error = ZSTD_isError(output_size);
        if (error > 0){
            std::cerr << "Decompression Error: - " << ZSTD_getErrorName(error) << std::endl;
            return -1;
        } */
        

        // free(compressed_buffer);
    }
}

/* namespace utility{
    std::string compress(char* data){
        namespace bio = boost::iostreams;
        
        std::stringstream compressed;
        std::stringstream origin;
        origin << data; //(data);

        bio::filtering_streambuf<bio::input> out;
        out.push(bio::zstd_compressor(bio::zstd_params(bio::zstd::default_compression)));

        out.push(origin);
        bio::copy(out, compressed);

        return compressed.str();
    }
    
    std::string decompress(char* data){
        namespace bio = boost::iostreams;
        
        std::stringstream decompressed;
        std::stringstream origin;
        origin << data; //(data);

        bio::filtering_streambuf<bio::input> out;
        out.push(bio::zstd_decompressor());

        out.push(origin);
        bio::copy(out, decompressed);

        return decompressed.str();
    } 
} */