#ifndef LAUNCHAPP_H
#define LAUNCHAPP_H

#include "launch.h"
#include <gtkmm.h>

class LaunchAppWindow;

class LaunchApplication: public Gtk::Application
{
    LaunchData &ld;
protected:
    LaunchApplication(LaunchData &ld);

public:
    static Glib::RefPtr<LaunchApplication> create(LaunchData &ld);

protected:
    void on_activate() override;

private:
    LaunchAppWindow* create_appwindow();
};

#endif
