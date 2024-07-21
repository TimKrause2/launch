CPPFLAGS=-ggdb -O2 -pthread -I/usr/include/PSOPT `pkg-config\
    --cflags freetype2 eigen3 ipopt gtkmm-4.0 epoxy` `sdl2-config --cflags`
LDLIBS=-L/usr/lib/PSOPT -lPSOPT -lGLESv2 -ljpeg -ltiff\
    `pkg-config --libs freetype2 ipopt adolc gtkmm-4.0 epoxy` `sdl2-config --libs`
LDFLAGS=-ggdb -pthread
CC=g++

COMMON_OBJS=gsim.o caams.o orbit.o image.o esfont.o shader.o uv_sphere.o esShader.o
SATELLITE_OBJS=satellite.o maneuver.o event.o sequencer.o two_phase_2.o thrust_vector.o\
    delayed_launch.o spheroid.o ballistic2.o icbm_simple.o

all:launch

orbit_test:orbit_test.o orbit.o caams.o

orbit_test.o:orbit_test.cpp

tracker:tracker.o $(COMMON_OBJS) tracker.o

launch:launch.o $(COMMON_OBJS) $(SATELLITE_OBJS) launchapp.o launchappwindow.o ballistic_dialog.o dialog_util.o

lunar:lunar.o $(COMMON_OBJS) $(SATELLITE_OBJS) trajectory.o


#	$(CC) $(LDFLAGS) lunar.o $(COMMON_OBJS) $(SATELLITE_OBJS) trajectory.o $(LDLIBS) -o lunar

orbitvis:orbitvis.o $(COMMON_OBJS)

gsim.o:gsim.cpp

caams.o:caams.cpp

orbit.o:orbit.cpp

image.o:image.c

font.o:font.cpp

esfont.o:esfont.cpp

esShader.o:esShader.c

sequencer.o:sequencer.cpp

orbitvis.o:orbitvis.cpp

launch.o:launch.cpp

lunar.o:lunar.cpp

trajectory.o:trajectory.cpp

shader.o:shader.cpp

uv_sphere.o:uv_sphere.cpp

satellite.o:satellite.cpp

maneuver.o:maneuver.cpp

event.o:event.cpp

two_phase.o:two_phase.cpp

two_phase_2.o:two_phase_2.cpp

thrust_vector.o:thrust_vector.cpp

tracker.o:tracker.cpp

delayed_launch.o:delayed_launch.cpp

spheroid.o:spheroid.cpp

ballistic.o:ballistic.cpp

ballistic2.o:ballistic2.cpp

icbm.o:icbm.cpp

icbm_simple.o:icbm_simple.cpp

icbm_simple2.o:icbm_simple2.cpp

icbm_simple_prototype.o:icbm_simple_prototype.cpp

launchapp.o:launchapp.cpp

launchappwindow.o:launchappwindow.cpp

ballistic_dialog.o:ballistic_dialog.cpp

dialog_util.o:dialog_util.cpp

clean:
	rm *.o
