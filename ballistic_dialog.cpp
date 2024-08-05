#include "ballistic_dialog.h"
#include "dialog_util.h"
#include <ios>
#include <iomanip>
#include <iostream>
#include <stdlib.h>
#include <cmath>


#define RAD_PER_DEG (M_PI/180.0)
#define DEG_PER_RAD (180.0/M_PI)
#define BL_LAT_VEHICLE  48.935333180089664*RAD_PER_DEG
#define BL_LONG_VEHICLE -54.58524740913247*RAD_PER_DEG
#define BL_LAT_TARGET   35.94248310361987*RAD_PER_DEG
#define BL_LONG_TARGET  -5.617076200801353*RAD_PER_DEG
#define BL_T1           500
#define BL_T2           1800
#define THRUST          (9.819609*2.0)

BallisticLaunchData::BallisticLaunchData()
{
    lat_launch = BL_LAT_VEHICLE;
    long_launch = BL_LONG_VEHICLE;
    lat_target = BL_LAT_TARGET;
    long_target = BL_LONG_TARGET;
    T1 = BL_T1;
    T2 = BL_T2;
    thrust = THRUST;
}

BallisticDialog::BallisticDialog(Gtk::Window& parent, BallisticLaunchData &bldata)
    :
      bldata(bldata),
      m_LabelLaunch("Launch site coordinates (latitude, longitude)"),
      m_LabelTarget("Target coordinates (latitude, longitude)"),
      m_LabelT1("Burn time (seconds)"),
      m_LabelT2("Coast time estimate (seconds)"),
      m_LabelThrust("Thrust acceleration (m/s/s)"),
      m_ButtonBox(Gtk::Orientation::HORIZONTAL, 0),
      m_Button_OK("Schedule"),
      m_Button_Cancel("Cancel"),
      gpsAlertMsg("Unrecognized GPS coordinate: "),
      numberAlertMsg("Unrecognized number: ")
{
    set_transient_for(parent);
    set_destroy_with_parent(true);

    auto provider = new_entry_provider();

    entry_set_font(m_EntryLaunch, provider);
    entry_set_font(m_EntryTarget, provider);
    entry_set_font(m_EntryT1, provider);
    entry_set_font(m_EntryT2, provider);
    entry_set_font(m_EntryThrust, provider);

    set_title("Ballistic Launch");
    set_child(m_Grid);
    m_Grid.set_expand(true);

    m_Grid.attach(m_LabelLaunch, 0, 0);
    m_Grid.attach(m_EntryLaunch, 1, 0);
    m_EntryLaunch.set_text(gps_to_text(bldata.lat_launch, bldata.long_launch));
    m_EntryLaunch.set_hexpand(true);

    m_Grid.attach(m_LabelTarget, 0, 1);
    m_Grid.attach(m_EntryTarget, 1, 1);
    m_EntryTarget.set_text(gps_to_text(bldata.lat_target, bldata.long_target));
    m_EntryTarget.set_hexpand(true);

    m_Grid.attach(m_LabelT1, 0, 2);
    m_Grid.attach(m_EntryT1, 1, 2);
    m_EntryT1.set_text(number_to_text(bldata.T1));
    m_EntryT1.set_hexpand(true);

    m_Grid.attach(m_LabelT2, 0, 3);
    m_Grid.attach(m_EntryT2, 1, 3);
    m_EntryT2.set_text(number_to_text(bldata.T2));
    m_EntryT2.set_hexpand(true);

    m_Grid.attach(m_LabelThrust, 0, 4);
    m_Grid.attach(m_EntryThrust, 1, 4);
    m_EntryThrust.set_text(number_to_text(bldata.thrust));
    m_EntryThrust.set_hexpand(true);

    m_Grid.attach(m_ButtonBox, 0, 5, 2, 1);
    m_ButtonBox.set_halign(Gtk::Align::END);
    m_ButtonBox.append(m_Button_OK);
    m_ButtonBox.append(m_Button_Cancel);

    m_EntryLaunch.signal_activate().connect(sigc::mem_fun(*this,&BallisticDialog::on_launch_activate));
    m_EntryTarget.signal_activate().connect(sigc::mem_fun(*this,&BallisticDialog::on_target_activate));
    m_EntryT1.signal_activate().connect(sigc::mem_fun(*this,&BallisticDialog::on_T1_activate));
    m_EntryT2.signal_activate().connect(sigc::mem_fun(*this,&BallisticDialog::on_T2_activate));
    m_EntryThrust.signal_activate().connect(sigc::mem_fun(*this,&BallisticDialog::on_thrust_activate));

}

BallisticDialog::~BallisticDialog()
{
}

void BallisticDialog::buttons_clicked_connect(
        const sigc::slot<void(const Glib::ustring&)>& slot)
{
    m_Button_OK.signal_clicked().connect(sigc::bind(slot, "OK"));
    m_Button_Cancel.signal_clicked().connect(sigc::bind(slot, "Cancel"));
}

void BallisticDialog::on_launch_activate()
{
    auto text = m_EntryLaunch.get_text();
    if(valid_gps(text)){
        text_to_gps(text, bldata.lat_launch, bldata.long_launch);
    }else{
        auto dialog = Gtk::AlertDialog::create(gpsAlertMsg + "Launch site");
        dialog->set_modal(true);
        dialog->show(*this);
        m_EntryLaunch.grab_focus();
    }
}

void BallisticDialog::on_target_activate()
{
    auto text = m_EntryTarget.get_text();
    if(valid_gps(text)){
        text_to_gps(text, bldata.lat_target, bldata.long_target);
    }else{
        auto dialog = Gtk::AlertDialog::create(gpsAlertMsg + "Target");
        dialog->set_modal(true);
        dialog->show(*this);
        m_EntryTarget.grab_focus();
    }
}

void BallisticDialog::on_T1_activate()
{
    auto text = m_EntryT1.get_text();
    if(valid_number(text)){
        text_to_number(text, bldata.T1);
    }else{
        auto dialog = Gtk::AlertDialog::create(numberAlertMsg + "T1");
        dialog->set_modal(true);
        dialog->show(*this);
        m_EntryT1.grab_focus();
    }
}

void BallisticDialog::on_T2_activate()
{
    auto text = m_EntryT2.get_text();
    if(valid_number(text)){
        text_to_number(text, bldata.T2);
    }else{
        auto dialog = Gtk::AlertDialog::create(numberAlertMsg + "T2");
        dialog->set_modal(true);
        dialog->show(*this);
        m_EntryT2.grab_focus();
    }
}

void BallisticDialog::on_thrust_activate()
{
    auto text = m_EntryThrust.get_text();
    if(valid_number(text)){
        text_to_number(text, bldata.thrust);
    }else{
        auto dialog = Gtk::AlertDialog::create(numberAlertMsg + "Thrust");
        dialog->set_modal(true);
        dialog->show(*this);
        m_EntryThrust.grab_focus();
    }
}
