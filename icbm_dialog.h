#include <gtkmm.h>

struct ICBMLaunchData
{
    double lat_launch, long_launch;
    double hgt_launch;
    double lat_target, long_target;
    double hgt_target;
    ICBMLaunchData();
};

class ICBMDialog : public Gtk::Window
{
private:
    ICBMLaunchData &icbm_data;
public:
    ICBMDialog(Gtk::Window& parent, ICBMLaunchData &icbm_data);
    ~ICBMDialog() override;

    void buttons_clicked_connect(const sigc::slot<void(const Glib::ustring&)>& slot);

protected:
    Gtk::Grid m_Grid;
    Gtk::Label m_LabelLaunch;
    Gtk::Label m_LabelLaunchHgt;
    Gtk::Label m_LabelTarget;
    Gtk::Label m_LabelTargetHgt;
    Gtk::Entry m_EntryLaunch;
    Gtk::Entry m_EntryLaunchHgt;
    Gtk::Entry m_EntryTarget;
    Gtk::Entry m_EntryTargetHgt;
    Gtk::Box m_ButtonBox;
    Gtk::Button m_Button_OK;
    Gtk::Button m_Button_Cancel;

    void on_launch_activate();
    void on_launch_hgt_activate();
    void on_target_activate();
    void on_target_hgt_activate();
};
