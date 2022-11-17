#include "scene.h"
#include <QMimeData>
#include <QGraphicsSceneDragDropEvent>
#include <QGraphicsSceneMouseEvent>
#include "resizableellipseitem.h"
#include "resizablepixmapitem.h"
#include "resizablerectitem.h"
#include "resizablestaritem.h"
#include <QDebug>

Scene::Scene(QObject *parent) : QGraphicsScene(parent),
    tool(Cursor), drawing(false), lastEraserCircle(nullptr)
{
    horGuideLine = addLine(-400,0,400,0,QPen(Qt::blue));
    verGuideLine = addLine(0,-400,0,400,QPen(Qt::blue));
    setSceneRect(-800,-400,1600,800);
}

void Scene::dragMoveEvent(QGraphicsSceneDragDropEvent *event)
{
    if(event->mimeData()->property("Key").canConvert(QMetaType::Int)){
        event->acceptProposedAction();
    }else{
        QGraphicsScene::dragMoveEvent(event);
    }
}

void Scene::dropEvent(QGraphicsSceneDragDropEvent *event)
{
    if(event->mimeData()->property("Key").canConvert(QMetaType::Int)){

        int key = event->mimeData()->property("Key").toInt();


        switch (key) {
        case 10:{
            //Ellipse
            ResizableEllipseItem * ellipse = new ResizableEllipseItem();
            ellipse->setRect(0,0,80,50);
            ellipse->setFlags(QGraphicsItem::ItemIsMovable | QGraphicsItem::ItemIsSelectable);
            ellipse->setBrush(Qt::gray);
            addItem(ellipse);

            ellipse->setPos(event->scenePos() -QPointF((ellipse->boundingRect().width()/2),
                                                       (ellipse->boundingRect().height()/2))) ;

        }
            break;
        case 20:{
            //Qt Quick Image
            ResizablePixmapItem * pixItem = new ResizablePixmapItem(QPixmap(":/IMAGES/wifi.png"));
            pixItem->setFlags(QGraphicsItem::ItemIsMovable | QGraphicsItem::ItemIsSelectable);
            addItem(pixItem);
            pixItem->setPos(event->scenePos() -QPointF((pixItem->boundingRect().width()/2),
                                                       (pixItem->boundingRect().height()/2))) ;
        }
            break;
        case 30:{
            //Rectangle
            ResizableRectItem * rectItem = new ResizableRectItem();
            rectItem->setRect(0,0,80,50);
            rectItem->setFlags(QGraphicsItem::ItemIsMovable | QGraphicsItem::ItemIsSelectable | QGraphicsItem::ItemIsFocusable);
            rectItem->setBrush(Qt::gray);
            addItem(rectItem);
            rectItem->setPos(event->scenePos() -QPointF((rectItem->boundingRect().width()/2),
                                                        (rectItem->boundingRect().height()/2))) ;
        }
            break;
        case 40:{
            //Star
            ResizableStarItem * starItem = new ResizableStarItem();
            starItem->setFlags(QGraphicsItem::ItemIsMovable | QGraphicsItem::ItemIsSelectable);
            starItem->setBrush(Qt::gray);
            addItem(starItem);
            starItem->setPos(event->scenePos() -QPointF((starItem->boundingRect().width()/2),
                                                        (starItem->boundingRect().height()/2))) ;
        }
            break;
        case 50:{
            //Qt Quick Image
            ResizablePixmapItem * pixItem1 = new ResizablePixmapItem(QPixmap(":/IMAGES/quick.png"));
            pixItem1->setFlags(QGraphicsItem::ItemIsMovable | QGraphicsItem::ItemIsSelectable);
            addItem(pixItem1);
            pixItem1->setPos(event->scenePos() -QPointF((pixItem1->boundingRect().width()/2),
                                                        (pixItem1->boundingRect().height()/2))) ;
        }
            break;

        }




        event->acceptProposedAction();
    }else{
        QGraphicsScene::dropEvent(event);
    }

}

void Scene::mousePressEvent(QGraphicsSceneMouseEvent *event)
{
    if(event->button() == Qt::LeftButton){
        if(tool == ToolType::Pen || tool == Eraser){
            //qDebug() << "Press event : " << event->button();
            startingPoint = event->scenePos();
            drawing = true;
        }else{
            QGraphicsScene::mousePressEvent(event);
        }
    }else{
        QGraphicsScene::mousePressEvent(event);
    }
}

void Scene::mouseMoveEvent(QGraphicsSceneMouseEvent *event)
{

    //if((event->button() == Qt::LeftButton) && drawing){ // not working
    if((event->buttons() & Qt::LeftButton) && drawing){
        //qDebug() << "Move event : " << event->button();
        if(tool == ToolType::Pen){
            //qDebug() << "Move event : " << event->button();
            drawLineTo(event->scenePos());
        }else if(tool == ToolType::Eraser){
            drawEraserAt(event->scenePos());
        }
    }else{
        QGraphicsScene::mouseMoveEvent(event);
    }
}

void Scene::mouseReleaseEvent(QGraphicsSceneMouseEvent *event)
{
    if((event->button() == Qt::LeftButton) && drawing){
        if(tool == ToolType::Pen){
            //qDebug() << "Release event : " << event->button();
            lineGroup = nullptr;
            drawing = false;
        }

        if(tool == ToolType::Eraser){
            removeItem(lastEraserCircle);
            delete lastEraserCircle;
            lastEraserCircle = nullptr;
            drawing = false;
        }


    }else{
        QGraphicsScene::mouseReleaseEvent(event);
    }
}

void Scene::drawLineTo(const QPointF &endPoint)
{
    if(!lineGroup){
        lineGroup = new QGraphicsItemGroup();
        lineGroup->setFlags(QGraphicsItem::ItemIsMovable | QGraphicsItem::ItemIsSelectable);
        addItem(lineGroup);
        lastPenPoint = startingPoint;
    }

    QGraphicsLineItem * localLine = new QGraphicsLineItem(QLineF(lastPenPoint, endPoint));
    QPen mPen;
    mPen.setWidth(3);
    mPen.setColor(Qt::green);
    localLine->setPen(mPen);
    lineGroup->addToGroup(localLine);

    lastPenPoint = endPoint;
}

void Scene::drawEraserAt(const QPointF &endPoint)
{
    if(!lastEraserCircle){
        lastEraserCircle = addEllipse(0,0,50,50);
    }
    lastEraserCircle->setPos(endPoint - QPointF(lastEraserCircle->boundingRect().width()/2,
                                                lastEraserCircle->boundingRect().height()/2));
    eraseStrokeUnder(lastEraserCircle);
}

void Scene::eraseStrokeUnder(QGraphicsEllipseItem *item)
{

    QList<QGraphicsItem *> itemsToRemove = item->collidingItems();
    QList<QGraphicsItemGroup *> groupItems;

    foreach(QGraphicsItem * myItem, itemsToRemove){

        QGraphicsItemGroup * group = dynamic_cast<QGraphicsItemGroup *>(myItem);
        if(group){
            groupItems.append(group);
        }

        //cast to graphicsLineItem
        QGraphicsLineItem * line = dynamic_cast<QGraphicsLineItem *>(myItem);
        if(line && (line != horGuideLine) && (line != verGuideLine)){
            qDebug() << "Group item has no child. Removing it";
            removeItem(line);
            delete line;
        }

    }

    //Remove group items that don't have any children.
    foreach(QGraphicsItemGroup * group, groupItems){
        if(group->childItems().count() == 0){
            removeItem(group);
            delete group;
        }
    }
}

Scene::ToolType Scene::getTool() const
{
    return tool;
}

void Scene::setTool(ToolType newTool)
{
    tool = newTool;
}
