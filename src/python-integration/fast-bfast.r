EnableFastBfast = function()
{
    if (exists("set_fast_options"))
    {
        print("Using fast BFAST")
        set_fast_options()
        # Disable bfastts modifications due to issue #2
        options(bfast.use_bfastts_modifications=FALSE)
    } else print("Using reference BFAST, install bfast2/bfast for a speed increase")
}
