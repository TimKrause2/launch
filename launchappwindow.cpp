#include "launchappwindow.h"
#include <iostream>

using std::cerr;
using std::endl;

LaunchAppWindow::LaunchAppWindow(LaunchData &ld)
    :
      Gtk::ApplicationWindow(),
      ld(ld),
      m_VBox(Gtk::Orientation::VERTICAL, 0),
      m_Button("test")
{
    set_title("Maneuver Tester");
    set_default_size(800,600);

    set_child(m_VBox);

    m_VBox.append(m_Button);

    m_GLArea.set_expand(true);
    m_GLArea.set_size_request(320,200);
    m_GLArea.set_auto_render(true);
    m_GLArea.set_use_es(true);
    m_GLArea.set_required_version(3,2);
    m_GLArea.set_has_depth_buffer(true);
    m_VBox.append(m_GLArea);

    m_GLArea.signal_realize().connect(sigc::mem_fun(*this,&LaunchAppWindow::gla_realize));
    m_GLArea.signal_unrealize().connect(sigc::mem_fun(*this,&LaunchAppWindow::gla_unrealize), false);
    m_GLArea.signal_render().connect(sigc::mem_fun(*this,&LaunchAppWindow::gla_render), false);
    m_GLArea.signal_resize().connect(sigc::mem_fun(*this,&LaunchAppWindow::gla_resize));

    auto drag = Gtk::GestureDrag::create();
    drag->signal_drag_begin().connect(sigc::mem_fun(*this,&LaunchAppWindow::drag_begin));
    drag->signal_drag_update().connect(sigc::mem_fun(*this,&LaunchAppWindow::drag_update));
    m_GLArea.add_controller(drag);

    auto escroll = Gtk::EventControllerScroll::create();
    escroll->set_flags(
                Gtk::EventControllerScroll::Flags::BOTH_AXES |
                Gtk::EventControllerScroll::Flags::DISCRETE);
    escroll->signal_scroll().connect(sigc::mem_fun(*this,&LaunchAppWindow::scroll),true);
    m_GLArea.add_controller(escroll);

    auto keys = Gtk::EventControllerKey::create();
    keys->signal_key_pressed().connect(sigc::mem_fun(*this,&LaunchAppWindow::key_pressed),false);
    add_controller(keys);

    m_timer_con = Glib::signal_timeout().connect(sigc::mem_fun(*this,&LaunchAppWindow::timeout), 17);

}

LaunchAppWindow::~LaunchAppWindow()
{
}

void LaunchAppWindow::gla_realize()
{
    m_GLArea.make_current();
    try
    {
        m_GLArea.throw_if_error();
        ld.init();
    }
    catch(const Gdk::GLError& gle)
    {
        cerr << "An error occured making the context current during realize:" << endl;
        cerr << gle.domain() << "-" << gle.code() << "-" << gle.what() << endl;
    }
}

void LaunchAppWindow::gla_unrealize()
{
    m_GLArea.make_current();
    try
    {
        m_GLArea.throw_if_error();
    }
    catch(const Gdk::GLError& gle)
    {
        cerr << "An error occured making the context current during unrealize" << endl;
        cerr << gle.domain() << "-" << gle.code() << "-" << gle.what() << endl;
    }
}

bool LaunchAppWindow::gla_render(const Glib::RefPtr<Gdk::GLContext>& /* context */)
{
    try
    {
        m_GLArea.throw_if_error();

        ld.display();

        glFlush();

        return true;
    }
    catch(const Gdk::GLError& gle)
    {
        cerr << "An error occurred in the render callback of the GLArea" << endl;
        cerr << gle.domain() << "-" << gle.code() << "-" << gle.what() << endl;
        return false;
    }
}

void LaunchAppWindow::gla_resize(int width, int height)
{
    ld.update_viewport(width, height);
}

void LaunchAppWindow::drag_begin(double start_x, double start_y)
{
    offset_x_last = 0.0;
    offset_y_last = 0.0;
}

void LaunchAppWindow::drag_update(double offset_x, double offset_y)
{
    double dx = offset_x - offset_x_last;
    double dy = offset_y - offset_y_last;
    ld.update_view_direction((int)dx, (int)dy);
    offset_x_last = offset_x;
    offset_y_last = offset_y;
}

bool LaunchAppWindow::scroll(double dx, double dy)
{
    ld.update_view_distance((int)dy);
    return true;
}

bool LaunchAppWindow::key_pressed(guint keyval, guint keycode, Gdk::ModifierType state)
{
    switch(keyval){
    case GDK_KEY_Up:
        ld.update_time_factor(1);
        return true;
    case GDK_KEY_Down:
        ld.update_time_factor(-1);
        return true;
    case GDK_KEY_F9:
        ballistic_dialog_show();
        return true;
    default:
        return false;
    }
}

bool LaunchAppWindow::timeout()
{
    m_GLArea.queue_draw();
    return true;
}

void LaunchAppWindow::ballistic_dialog_show()
{
    BallisticDialog* pDialog = new BallisticDialog(*this, bldata);
    pDialog->set_modal(true);
    pDialog->buttons_clicked_connect(
                sigc::bind(
                    sigc::mem_fun(*this, &LaunchAppWindow::on_ballistic_dialog_response), pDialog));
    pDialog->signal_close_request().connect(sigc::bind(
        sigc::mem_fun(*this, &LaunchAppWindow::on_ballistic_dialog_close), pDialog), false);
    pDialog->set_visible(true);
}

void LaunchAppWindow::on_ballistic_dialog_response(
        const Glib::ustring& response,
        BallisticDialog* dialog)
{
    delete dialog;

    if(response == "OK"){
        if(!ld.body_sat->ScheduleICBMLaunch(
                    bldata.lat_launch, bldata.long_launch,
                    bldata.lat_target, bldata.long_target,
                    bldata.T1, bldata.T2,
                    bldata.thrust))
        {
            auto dialog = Gtk::AlertDialog::create("Optimizer failed!!! See terminal output.");
            dialog->set_modal(true);
            dialog->show(*this);
        }
    }
}

bool LaunchAppWindow::on_ballistic_dialog_close(BallisticDialog *pDialog)
{
    delete pDialog;
    return false;
}
