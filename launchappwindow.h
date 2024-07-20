#ifndef LAUNCHAPPWINDOW_H
#define LAUNCHAPPWINDOW_H

#include "launch.h"
#include "ballistic_dialog.h"
#include <gtkmm.h>

class LaunchAppWindow : public Gtk::ApplicationWindow
{
private:
    LaunchData &ld;
    BallisticLaunchData bldata;
    double offset_x_last;
    double offset_y_last;
public:
    LaunchAppWindow(LaunchData &ld);
    ~LaunchAppWindow() override;
protected:
    Gtk::Box m_VBox;
    Gtk::Button m_Button;
    Gtk::GLArea m_GLArea;
    sigc::connection m_timer_con;

    void gla_realize();
    void gla_unrealize();
    bool gla_render(const Glib::RefPtr<Gdk::GLContext>& context);
    void gla_resize(int width, int height);
    void drag_begin(double start_x, double start_y);
    void drag_update(double offset_x, double offset_y);
    bool scroll(double dx, double dy);
    bool key_pressed(guint keyval, guint keycode, Gdk::ModifierType state);
    bool timeout();
    void ballistic_dialog_show();
    void on_ballistic_dialog_response(const Glib::ustring& response, BallisticDialog* dialog);
    bool on_ballistic_dialog_close(BallisticDialog* dialog);
};

#endif
