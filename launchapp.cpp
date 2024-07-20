#include "launchapp.h"
#include "launchappwindow.h"

LaunchApplication::LaunchApplication(LaunchData &ld)
    :Gtk::Application("krause.tim.launch", Gio::Application::Flags::NON_UNIQUE),
      ld(ld)
{
}

Glib::RefPtr<LaunchApplication> LaunchApplication::create(LaunchData &ld)
{
    return Glib::make_refptr_for_instance<LaunchApplication>(new LaunchApplication(ld));
}

LaunchAppWindow *LaunchApplication::create_appwindow()
{
    auto appwindow = new LaunchAppWindow(ld);

    add_window(*appwindow);

    appwindow->signal_hide().connect([appwindow](){ delete appwindow; });

    return appwindow;
}

void LaunchApplication::on_activate()
{
    auto appwindow = create_appwindow();
    appwindow->present();
}



