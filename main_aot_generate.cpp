/*
 * 2014-12-18, 01:37
 * algorithm part is finished roughly.
 * Note that some refinements can be explored. Please search "TODO" labels.
 * 
 * Schdule part finished. It runs much faster than default schedule.
 * But still lots of things to try. Furthermore, this implementation is
 * much faster than Matlab implementation though Matlab implementation is 
 * complete haze-removal algorithm. Try to find better scheduling.
 *
 * Personal Notes:
 * on image size=[2488, 2488, 3]
 *
 *   MATLAB took 1.7xxx sec.
 *   This program took 0.155x sec if ALL compute_root
 *   This program took 0.283x sec if scheduling of finding parameters
 *   (I have comment them out for testing ALL compute_root)
 *
 * future work:
 * Refine algorithm to obtain complete implementation
 * Use AOT compiling for this program and porting to Android.
 *
 */

// Those two parameters have critical imapct on performance
#define STRIP_WIDTH 4 // width of parallel computing strip
#define VEC_LEN 16 // in x86 structure(SSE instructions), may use length 4*N?

#define ITER 5 // iterating times for measuring timing

#include <stdio.h>
#include <Halide.h>

using Halide::Image;
#include "image_io.h"

#include "clock.h"

using namespace Halide;

int main(int argc, char **argv){
    
    // color image
    ImageParam input(type_of<uint8_t>(), 3);
    //ImageParam input(type_of<float>(), 3);

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
    sigma2_n(0) = sum(sigma2_wnr(0+box_wimg.x, 0+box_wimg.y))/(input.width()*input.height());
    
    Func weight;
    weight(x, y) = max((sigma2_wnr(x, y) - sigma2_n(0)), 0.0f)/max(sigma2_n(0), sigma2_wnr(x, y));
    Func dcpf_refine;
    dcpf_refine(x, y) = mu_wnr(x, y) + weight(x, y)*(dcpf_p(x, y)-mu_wnr(x, y));

    // get transmission map
    Func t_map;
    t_map(x, y) = clamp(1.0f-0.95f*dcpf_refine(x, y), 0.01f, 1.0f);

    // recover image
    Func restore_img;
    restore_img(x, y, c) = cast<uint8_t>(clamp((clamped(x, y, c)-airlight(c))/t_map(x, y) + airlight(c), 0.0f, 1.0f)*255.0f);


    /*
     * Schedule Parts 
     */
    // TODO: need to experiment various schedule for haze removal

    //    
    // begin to filter DCP with adaptive wiener filter
    // from bottom(restore_img) to up(dcp_f) schduling
    //

    // Now we can sigma2_n.compute_root() and do parallel both at restore_img and inside sigma2_n
    // It's becuase calculating sigma2_n has RDom box_wimg. So to avoid to scan whole
    // image for every parallel for-loop, we need sigma2_wnr and sigma2_n stored at root.
    //
    // Strategy one: divide algorithm to 2 parts:
    // Part1 for getting all wienre filter paramters and store at root. Parallel can be applied when 
    // calculating paramters.
    // Part2 for restoring image. Parallel can be applied, too.

    /* restore_img */
    // compute color channels at innermost for-loopt
    restore_img.reorder(c, x, y)
        .bound(c, 0, 3)
        .unroll(c);

    // Parallel here require at least sigma2_n stored at root.
    // Or we may need to scan whole image for every parallel for-loop because sigma2_n need whole simga2_wnr.
    // to be faster, it'll be ok to store mu_wnr, sigma2_wnr at root. But for saving memory, it's also ok to compute 
    // sigma2_wnr and mu_wnr here again. Note that we have already calculated sigma2_wnr and mu_wnr one time for getting
    // sigma2_n

    restore_img.parallel(y);
    restore_img.vectorize(x, VEC_LEN);

    /* clapmed, airlight, t_map */
    // No neighbor computing in path t_map->dcpf_refine->weight, so all inlined into restore_img

    // clamped is used almost all stages. Even airlight needs WHOLE clamped, so compute_root()
    clamped.reorder(c, x ,y)
        .bound(c, 0, 3)
        .unroll(c);
    clamped.parallel(y);
    clamped.vectorize(x, VEC_LEN);
    clamped.compute_root();
    // TODO: test no reorder 'c' version. Then we can use airlight.compute_at(restore_img, c)
    airlight.parallel(c);
    airlight.compute_root();
   
    /* dcpf_p */
    dcpf_p.parallel(y);
    dcpf_p.vectorize(x, VEC_LEN);
    dcpf_p.compute_root();
    
    /* clamped_f_no_air */
    // this is only used in calculating dcpf_p. No need to store
    // and keep inline

    //
    // begin to extract wiener filter paramters. 
    // from bottom(sigma2_n) to up(local_sum_dcp and local_sum2_dcp)
    //

    Var yi, yo; // for neighboring computing
    /* sigma2_n */
    sigma2_n.compute_root();

    /* sigma2_wnr */
    /*
    sigma2_wnr.split(y, yo, yi, STRIP_WIDTH).parallel(yo);
    sigma2_wnr.vectorize(x, VEC_LEN);
    sigma2_wnr.compute_root();
    */
    // test all compute_root()
    sigma2_wnr.compute_root().parallel(y).vectorize(x, VEC_LEN);

    /* local_sum2_dcp */
    /*
    local_sum2_dcp.store_at(sigma2_wnr, yo).compute_at(sigma2_wnr, yi);
    local_sum2_dcp.vectorize(x, VEC_LEN);
    local_sum2_dcp.update(0).split(y, yo, yi, STRIP_WIDTH);
    local_sum2_dcp.update(0).vectorize(x, VEC_LEN);
    */
    // test all compute_root()
    local_sum2_dcp.compute_root().parallel(y).vectorize(x, VEC_LEN);
    local_sum2_dcp.update(0).parallel(y).vectorize(x, VEC_LEN);
    /* dcpf2 */
    /*
    dcpf2.store_at(local_sum2_dcp, yo).compute_at(local_sum2_dcp, yi);
    dcpf2.vectorize(x, VEC_LEN);
    */
    // test all compute_root()
    dcpf2.compute_root().parallel(y).vectorize(x, VEC_LEN);

    /* mu_wnr */
    /*
    mu_wnr.split(y, yo, yi, STRIP_WIDTH).parallel(yo);
    mu_wnr.vectorize(x, VEC_LEN);
    mu_wnr.compute_root();
    */
    // test all compute_root()
    mu_wnr.compute_root().parallel(y).vectorize(x, VEC_LEN);
    /* local_sum_dcp */
    /*
    local_sum_dcp.store_at(mu_wnr, yo).compute_at(mu_wnr, yi);
    local_sum_dcp.vectorize(x, VEC_LEN);
    local_sum_dcp.update(0).vectorize(x, VEC_LEN);
    */
    // test all compute_root()
    local_sum_dcp.compute_root().parallel(y).vectorize(x, VEC_LEN);
    local_sum_dcp.update(0).parallel(y).vectorize(x, VEC_LEN);

    /*
     * AOT generate
     */

    std::vector<Argument> args(1);
    args[0] = input;
    restore_img.compile_to_file("halide_haze_removal", args);

    printf("Halide pipeline generated, but not yet run.\n");

    return 0;
}// int main
