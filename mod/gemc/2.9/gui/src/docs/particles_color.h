QString PColorsTab::load_doc()
{
	return
	"<html><body>"
	"<br>"
	"<center>"
	"<h2> Paricle Color Codes</h2>"
	"</center>"
	"<br><br>"
	"<table bordersize='0'>"
	"<tr><td>Particle</td>            <td>Color</td></tr>"
	"<tr><td><br></td><td></td></tr>"
	"<tr><td>neutron:</td>           <td><font color='black'>  black </font></td></tr>"
	"<tr><td>gamma:</td>            <td><font color='blue'>   blue   </font></td></tr>"
        "<tr><td>optical photon:</td>            <td><font color='black'>   white   </font></td></tr>"        
	"<tr><td>e-:  </td>               <td><font color='red'>   red   </font></td></tr>"
        "<tr><td>e+:  </td>               <td><font color='yellow'>   yellow   </font></td></tr>"        
	"<tr><td>proton: </td>           <td><font color='orange'> orange </font></td></tr>"
	"<tr><td>&pi;<sup>+</sup>/&pi;<sup>-</sup>:</td>   <td><font color='magenta'> magenta </font></td></tr>"
	"<tr><td>&mu;<sup>-</sup>/&mu;<sup>-</sup>:</td>   <td><font color='brown'>  brown  </font></td></tr>"
	"<tr><td><br></td><td></td></tr>"
	"<tr><td colspan='2'>all other particle colors are based on their charges:</td></tr>"
	"<tr><td><br></td><td></td></tr>"
	"<tr><td>q=+1:</td>               <td><font color='cyan'>   cyan   </font></td></tr>"
	"<tr><td>q= 0:</td>               <td><font color='gray'> gray </font></td></tr>"
	"<tr><td>q=-1:</td>               <td><font color='green'> green </font></td></tr>"
	"</table>"
	" </body></html>";
}
