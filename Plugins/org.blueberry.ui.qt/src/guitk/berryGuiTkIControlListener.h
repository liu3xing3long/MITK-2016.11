/*===================================================================

BlueBerry Platform

Copyright (c) German Cancer Research Center,
Division of Medical and Biological Informatics.
All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef BERRYGUITKICONTROLLISTENER_H_
#define BERRYGUITKICONTROLLISTENER_H_

#include <berryMacros.h>
#include <berryMessage.h>

#include <org_blueberry_ui_qt_Export.h>
#include "berryGuiTkControlEvent.h"

namespace berry
{

namespace GuiTk
{

/**
 * Classes which implement this interface provide methods
 * that deal with the events that are generated by moving
 * and resizing controls.
 * <p>
 * After creating an instance of a class that implements
 * this interface it can be added to a control using the
 * <code>addControlListener</code> method and removed using
 * the <code>removeControlListener</code> method. When a
 * control is moved or resized, the appropriate method will
 * be invoked.
 * </p>
 *
 * @see ControlAdapter
 * @see ControlEvent
 */
struct BERRY_UI_QT IControlListener: public virtual Object
{

  berryObjectMacro(berry::GuiTk::IControlListener)

  struct BERRY_UI_QT Events {

    enum Type {
     NONE      = 0x00000000,
     MOVED     = 0x00000001,
     RESIZED   = 0x00000002,
     ACTIVATED = 0x00000004,
     DESTROYED = 0x00000008,

     ALL       = 0xffffffff
    };

    Q_DECLARE_FLAGS(Types, Type)

    typedef Message1<ControlEvent::Pointer> EventType;

    EventType movedEvent;
    EventType resizedEvent;
    EventType activatedEvent;
    EventType destroyedEvent;

    void AddListener(IControlListener::Pointer listener);
    void RemoveListener(IControlListener::Pointer listener);

  private:
    typedef MessageDelegate1<IControlListener, ControlEvent::Pointer> Delegate;
  };

  virtual ~IControlListener();

  virtual Events::Types GetEventTypes() const = 0;

  /**
   * Sent when the location (x, y) of a control changes relative
   * to its parent (or relative to the display, for <code>Shell</code>s).
   *
   * @param e an event containing information about the move
   */
  virtual void ControlMoved(ControlEvent::Pointer /*e*/)
  {
  }

  /**
   * Sent when the size (width, height) of a control changes.
   *
   * @param e an event containing information about the resize
   */
  virtual void ControlResized(ControlEvent::Pointer /*e*/)
  {
  }

  virtual void ControlActivated(ControlEvent::Pointer /*e*/)
  {
  }

  virtual void ControlDestroyed(ControlEvent::Pointer /*e*/)
  {
  }

};

template<typename R>
struct ControlMovedAdapter: public IControlListener
{
  typedef R Listener;
  typedef void
      (R::*Callback)(ControlEvent::Pointer);

  ControlMovedAdapter(R* l, Callback c) :
    listener(l), callback(c)
  {
    poco_assert(listener);
    poco_assert(callback);
  }

  Events::Types GetEventTypes() const override
  {
    return Events::MOVED;
  }

  void ControlMoved(ControlEvent::Pointer e) override
  {
    (listener->*callback)(e);
  }

private:

  Listener* listener;
  Callback callback;
};

template<typename R>
struct ControlResizedAdapter: public IControlListener
{
  typedef R Listener;
  typedef void
      (R::*Callback)(ControlEvent::Pointer);

  ControlResizedAdapter(R* l, Callback c) :
    listener(l), callback(c)
  {
    poco_assert(listener);
    poco_assert(callback);
  }

  Events::Types GetEventTypes() const override
  {
    return Events::RESIZED;
  }

  void ControlResized(ControlEvent::Pointer e) override
  {
    (listener->*callback)(e);
  }

private:

  Listener* listener;
  Callback callback;
};

template<typename R>
struct ControlActivatedAdapter: public IControlListener
{
  typedef R Listener;
  typedef void
      (R::*Callback)(ControlEvent::Pointer);

  ControlActivatedAdapter(R* l, Callback c) :
    listener(l), callback(c)
  {
    poco_assert(listener);
    poco_assert(callback);
  }

  Events::Types GetEventTypes() const override
  {
    return Events::ACTIVATED;
  }

  void ControlActivated(ControlEvent::Pointer e) override
  {
    (listener->*callback)(e);
  }

private:

  Listener* listener;
  Callback callback;
};

template<typename R>
struct ControlDestroyedAdapter: public IControlListener
{
  typedef R Listener;
  typedef void
      (R::*Callback)(ControlEvent::Pointer);

  ControlDestroyedAdapter(R* l, Callback c) :
    listener(l), callback(c)
  {
    poco_assert(listener);
    poco_assert(callback);
  }

  Events::Types GetEventTypes() const override
  {
    return Events::DESTROYED;
  }

  void ControlDestroyed(ControlEvent::Pointer e) override
  {
    (listener->*callback)(e);
  }

private:

  Listener* listener;
  Callback callback;
};

}

}

Q_DECLARE_OPERATORS_FOR_FLAGS(berry::GuiTk::IControlListener::Events::Types)

#endif /* BERRYGUITKICONTROLLISTENER_H_ */
