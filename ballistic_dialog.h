#include <gtkmm.h>

struct BallisticLaunchData
{
    double lat_launch, long_launch;
    double lat_target, long_target;
    double T1, T2;
    double thrust;
    BallisticLaunchData();
};

class BallisticDialog : public Gtk::Window
{
private:
    BallisticLaunchData &bldata;
    Glib::ustring gpsAlertMsg;
    Glib::ustring numberAlertMsg;
public:
    BallisticDialog(Gtk::Window& parent, BallisticLaunchData &bldata);
    ~BallisticDialog() override;

    void buttons_clicked_connect(const sigc::slot<void(const Glib::ustring&)>& slot);

protected:
    Gtk::Grid m_Grid;
    Gtk::Label m_LabelLaunch;
    Gtk::Label m_LabelTarget;
    Gtk::Label m_LabelT1;
    Gtk::Label m_LabelT2;
    Gtk::Label m_LabelThrust;
    Gtk::Entry m_EntryLaunch;
    Gtk::Entry m_EntryTarget;
    Gtk::Entry m_EntryT1;
    Gtk::Entry m_EntryT2;
    Gtk::Entry m_EntryThrust;
    Gtk::Box m_ButtonBox;
    Gtk::Button m_Button_OK;
    Gtk::Button m_Button_Cancel;

    void on_launch_activate();
    void on_target_activate();
    void on_T1_activate();
    void on_T2_activate();
    void on_thrust_activate();
};
