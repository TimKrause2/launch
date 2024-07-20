#include "event.h"
#include "satellite.h"
#include <iostream>

Event::Event(TimeSpec t_absolute)
{
    t = t_absolute;
}

Event::~Event(){}

bool Event::Expired(TimeSpec &tnow)
{
    return t<=tnow;
}

SatelliteEvent::SatelliteEvent(Satellite* satellite, TimeSpec t_absolute)
    :Event(t_absolute),
      satellite(satellite)
{}

SatelliteEvent::~SatelliteEvent()
{}

IntegratorEvent::IntegratorEvent(TimeSpec t_absolute)
    : Event(t_absolute)
{}

void SatelliteEvent::Action(void)
{
    satellite->maneuverQueue.DeleteCurrent();
}

void EventQueue::AddEvent(Event *event)
{
    //std::cout << "AddEvent t = " << event->t.to_str() << std::endl;
    events.push_back(event);
}

bool EventQueue::Empty(void)
{
    return events.empty();
}

void EventQueue::Check(void)
{
    while(true)
    {
        if(events.empty()) return;
        Event* event = events.front();
        if(event->Expired(tnow)){
            events.pop_front();
            event->Action();
            //std::cout << "Deleting event. t = " << event->t.to_str() << std::endl;
            delete event;
        }else{
            break;
        }
    }
}

void EventQueue::Update(double dt)
{
    tnow += dt;
}

double EventQueue::NextEvent(void)
{
    return (events.front()->t - tnow).to_double();
}

TimeSpec EventQueue::AbsoluteTime(double time_from_now)
{
    return tnow + time_from_now;
}
