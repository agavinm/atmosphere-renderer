# Atmosphere renderer  
  
Atmosphere renderer developed in C++ using the [Mitsuba 2 renderer](https://github.com/mitsuba-renderer/mitsuba2).  

This Bachelor's Degree final project is based on [Guimera's physically-based optical model of the atmosphere](http://giga.cps.unizar.es/~ajarabo/pubs/Guimera2018spatio/index.html).

### Documentation
Report: [PDF](https://github.com/agavinm/atmosphere-renderer/files/5781595/Report.pdf) (in Spanish)  
Slides: [ODP](https://github.com/agavinm/atmosphere-renderer/files/5781603/Slides.zip) | [PDF](https://github.com/agavinm/atmosphere-renderer/files/5781604/Slides.pdf)  

### Renders
Render from outer space without aerosols (Pinhole camera). [XML](https://raw.githubusercontent.com/agavinm/atmosphere-renderer/master/scenes/space/space2.xml)  
<img src="https://user-images.githubusercontent.com/37375662/103893781-2221a180-50ee-11eb-8bb9-2bb3d88cd68b.png" width="640">  

Render from inside Earth's atmosphere without aerosols (Pinhole camera). The Sun is 90 degrees above the horizon. [XML](https://github.com/agavinm/atmosphere-renderer/files/5783353/render8.xml.txt)  
<img src="https://user-images.githubusercontent.com/37375662/103893786-264dbf00-50ee-11eb-8546-8d93c0469cc0.png" width="640">  

Render from inside Earth's atmosphere without aerosols (Pinhole camera). The Sun is 180 degrees above the horizon. [XML](https://github.com/agavinm/atmosphere-renderer/files/5783355/render17.xml.txt)  
<img src="https://user-images.githubusercontent.com/37375662/103893791-28b01900-50ee-11eb-98fc-de4837d0d3aa.png" width="640">  

Renders from inside Earth's atmosphere with Desert Dust environment (Fisheye camera) and different atmosphere turbidity (increases from left to right). The Sun is 90 degrees above the horizon. [XML](https://raw.githubusercontent.com/agavinm/atmosphere-renderer/master/scenes/latest/set1/renderF12A.xml) | [XML](https://raw.githubusercontent.com/agavinm/atmosphere-renderer/master/scenes/latest/set1/renderF12B.xml) | [XML](https://raw.githubusercontent.com/agavinm/atmosphere-renderer/master/scenes/latest/set1/renderF12C.xml)  
<img src="https://user-images.githubusercontent.com/37375662/103893821-36659e80-50ee-11eb-8a88-8c75f261fa74.png" width="160">
<img src="https://user-images.githubusercontent.com/37375662/103893823-382f6200-50ee-11eb-9496-fa90d41ed78f.png" width="160">
<img src="https://user-images.githubusercontent.com/37375662/103893830-39f92580-50ee-11eb-9ece-00133cfe6986.png" width="160">  

Renders from inside Earth's atmosphere with different aerosol models (Fisheye camera). The Sun is 30 degrees above the horizon. [XML](https://raw.githubusercontent.com/agavinm/atmosphere-renderer/master/scenes/latest/set4/renderF12C.xml) | [XML](https://raw.githubusercontent.com/agavinm/atmosphere-renderer/master/scenes/latest/set4/renderF12Ca.xml) | [XML](https://raw.githubusercontent.com/agavinm/atmosphere-renderer/master/scenes/latest/set4/renderF14C.xml) | [XML](https://raw.githubusercontent.com/agavinm/atmosphere-renderer/master/scenes/latest/set4/renderF16C.xml)  
<img src="https://user-images.githubusercontent.com/37375662/103893840-3fef0680-50ee-11eb-9904-33d87fda0096.png" width="160">
<img src="https://user-images.githubusercontent.com/37375662/103893847-41b8ca00-50ee-11eb-92c3-04a989afd3c4.png" width="160">
<img src="https://user-images.githubusercontent.com/37375662/103893855-44b3ba80-50ee-11eb-8d62-8c75fe7cb605.png" width="160">
<img src="https://user-images.githubusercontent.com/37375662/103893859-47aeab00-50ee-11eb-953b-57c418f116c9.png" width="160">

### About
* Author: Pedro Andrés Gavín Murillo  
* Supervisor: Adrián Jarabo Torrijos  
