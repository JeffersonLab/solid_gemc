// The permanent magnet is made by ring with cad file "../geometry/magnet_cad_11rings/MagnetCubes_21.stl" and ../geometry/magnet_cad_11rings/AnnularPlateHardcoded-AnnularPlate_21.stl"
// each ring is about 24mm long with 15mm diameter inner space
//its map is as follows (1st column is radius relative to main axis z, 2nd column is Bx field, while By=0 and Bz=0)
// 0.0000		1961.3753
// 3.7331		1964.1460
// 7.4662		1972.5013
// 11.1993		1986.5709
// 14.9324		2006.5761
// 18.6655		2032.8355
// 22.3986		2065.7753
// 26.1317		2105.9417
// 29.8648		2154.0180
// 33.5979		2210.8442
// 37.3311		2277.4421
// 41.0642		2355.0427
// 44.7973		2445.1179
// 48.5304		2549.4099
// 52.2635		2669.9555
// 55.9966		2809.0878
// 59.7297		2969.3885
// 63.4628		3153.5288
// 67.1959		3363.8778
// 70.9290		3601.6357

 // refer to "map_shape_20221022_1.png" and "map_shape_20221022_2.png"

{
ofstream output("output",ios::trunc);

    // this function can fit the data, refer to "map_fit_20251022.png"
    TF1 *f=new TF1("f","4.2659e-3 * pow(x,3) - 2.8327e-2 * pow(x,2) + 3.6691 * x + 1950.7985",0,70);

    for (double x=-70;x<=70;x=x+7){
        for (double y=-70;y<=70;y=y+7){
            // for (double z=0;z<=24;z=z+24){ //1 ring
            // for (double z=0;z<=254;z=z+254){ //11 rings
            // for (double z=0;z<=381;z=z+381){//16 rings
            for (double z=0;z<=813;z=z+813){//33 rings
                output << x << "\t" << y << "\t" << z << "\t" << f->Eval(sqrt(x*x+y*y)) << "\t" << 0 << "\t" << 0 << endl;
            }
        }
    }

output.close();
}
