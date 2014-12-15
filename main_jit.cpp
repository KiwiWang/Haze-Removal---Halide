/*
 * 2014-12-15, 20:45
 * algorithm part is finished roughly.
 * Note that some refinements can be explored. Please search "TODO" labels.
 *
 * However, this implementation is super slow(even slower than Matlab version)
 * It's because shedule part is not explicitly set. Now just use default
 * schedule, that is, 'inline' for almost all Func.
 *
 * So even though algorithm part needs some refinements, it may be more 
 * important to come up with a reasonable schedule part. And then we can 
 * use AOT compiing for this program or porting to Android.
 *
 */


#include <stdio.h>
#include <Halide.h>

using Halide::Image;
#include "image_io.h"

using namespace Halide;

int main(int argc, char **argv){
    
    const char *output_name = "output.png";
    if(argc==2){
        fprintf(stdout, "Output filename:%s\n", output_name);
    } else if(argc==3){
        fprintf(stdout, "Output filename:%s\n", argv[2]);
    } else{
        fprintf(stderr, "Usage: ./main_jit <input_file> [output_file]\n");
        return -1;
    }//if(argc==2)

    Image<uint8_t> input = load<uint8_t>(argv[1]);

    Var x, y, c;

    // set boundary condition
    Func clamped;
    Expr x_clamped = clamp(x, 0, input.width()-1);
    Expr y_clamped = clamp(y, 0, input.height()-1);

    // sqrt version does well with preserving color.
    // remove sqrt if so desire.
    clamped(x, y, c) = cast<float>(sqrt(input(x_clamped, y_clamped, c)/255.0f));
    //clamped(x, y, c) = cast<float>(input(x_clamped, y_clamped, c)/255.0f);
    
    /*
     *  algorithm parts  
     */ 

    // find air light
    // 1. here will do better if use argmax to find most heavily hazy location
    //    this appoach will need find DCP to detect heaviliness of haze first.
    //
    // 2. Or if we can sort whole image, select top 10% brightest pixels and 
    //    calculate their mean to get a better airlight estimation.
    //
    // However, now we just use simplest way to estimate airlight
    // from brightest pixel for each channel

    /** TODO:need to use argmax to find brightest location... **/

    // use rough conversion for Y channel
    //Func lut; // use look up table for gamma correction
    //Var i;
    //lut(i) = clamp(pow(i/255.0f, 2.2f)*255.0f, 0.0f, 255.5f);

    //Func gray_img;
    //gray_img(x, y) = 0.299f*lut(clamped(x, y, 0)) + 0.587f*lut(clamped(x, y, 1)) + 0.114f*lut(clamped(x, y, 2));
    //gray_img(x, y) = clamped(x, y, 0)/3.0f + clamped(x, y, 1)/3.0f + clamped(x, y, 2)/3.0f;
    

    // find brightest pixel for each channel
    RDom box_wimg(0, input.width()-1, 0, input.height()-1);
    // TODO: this is a simplest but very rough estimation for airligh.
    // fix it using argmax as mentioned above
    Func airlight;
    airlight(c) = maximum(clamped(0+box_wimg.x, 0+box_wimg.y, c));

    // filter out airlight
    Func clamped_f_no_air, dcpf_p;
    clamped_f_no_air(x, y, c) = clamped(x, y, c)/airlight(c);
    // and find DCP we want
    dcpf_p(x, y) = min(clamped_f_no_air(x, y, 0), min(clamped_f_no_air(x, y, 1), clamped_f_no_air(x, y, 2)));

    // begin to estimate Adaptive Wiener filter parameters
    // TODO: actually, this filter step should be applied twice.
    // Re-write this step as function. But need to think about schedule
    int kernel_size = 21; // may have better way to determine, but roughly 21
    RDom box_wnr(-kernel_size/2, kernel_size, -kernel_size/2, kernel_size);

    // TODO: is there more compact implementation? it's too many Func I think...
    
    Func local_sum_dcp, mu_wnr; // for finding mu
    local_sum_dcp(x, y) = 0.0f;
    local_sum_dcp(x, y) += dcpf_p(x + box_wnr.x, y + box_wnr.y);
    mu_wnr(x, y) = local_sum_dcp(x, y)/(float)(kernel_size*kernel_size);

    Func dcpf2, local_sum2_dcp, sigma2_wnr; // for calculating sigma2
    dcpf2(x, y) = dcpf_p(x, y)*dcpf_p(x, y);
    local_sum2_dcp(x, y) = 0.0f;
    local_sum2_dcp(x, y) += dcpf2(x + box_wnr.x, y + box_wnr.y);
    sigma2_wnr(x, y) = max(local_sum2_dcp(x, y)/(float)(kernel_size*kernel_size) - mu_wnr(x, y)*mu_wnr(x, y), 0.0f);

    // begin to filter DCP (dcpf_p) with adaptive wiener filter
    //TODO: is this a practical implementation to find sigma2_n?
    Func sigma2_n;
    Var dummy;
    sigma2_n(dummy) = 0.0f;
    sigma2_n(0) = sum(sigma2_wnr(0+box_wimg.x, 0+box_wimg.y))/(float)(input.width()*input.height());
    
    Func weight;
    weight(x, y) = max((sigma2_wnr(x, y) - sigma2_n(0)), 0.0f)/max(sigma2_n(0), sigma2_wnr(x, y));
    Func dcpf_refine;
    dcpf_refine(x, y) = mu_wnr(x, y) + weight(x, y)*(dcpf_p(x, y)-mu_wnr(x, y));

    // get transmission map
    Func t_map;
    t_map(x, y) = clamp(1.0f-0.95f*dcpf_refine(x, y), 0.01f, 1.0f);

    // recover image
    Func restore_img;
    restore_img(x, y, c) = clamp((clamped(x, y, c)-airlight(c))/t_map(x, y) + airlight(c), 0.0f, 1.0f);


    /*
     * Schedule Parts 
     */
    // TODO: need to experiment various schedule for haze removal
    // There are many stages and sub-stages in this implementation.
    // We may need to plot some diagrams to help finding better schduling.

    //lut.compute_root();
    airlight.compute_root();
    sigma2_n.compute_root();


    /*
     * JIT output
     */
    // Func to testing output...
    Func test;
    test(x, y, c) = cast<uint8_t>(restore_img(x, y, c)*255.0f);

    Image<uint8_t> output = test.realize(input.width(), input.height(), input.channels());
    if(argc==2){
        save(output, output_name);
        return 0;
    }else
        save(output, argv[2]);
    return 0;
}// int main
